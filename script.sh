#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Check if Rust is installed
if ! command_exists cargo; then
    echo "Rust is not installed. Please install Rust and Cargo first."
    echo "Visit https://www.rust-lang.org/tools/install for installation instructions."
    exit 1
fi

# Create a new Rust project
PROJECT_NAME="jsonl_processor"
cargo new "$PROJECT_NAME"
cd "$PROJECT_NAME"

# Update Cargo.toml with dependencies
cat << EOF >> Cargo.toml

serde = { version = "1.0", features = ["derive"] }
serde_json = "1.0"
flate2 = "1.0"
EOF

# Create the main.rs file with our updated JSONL processor code
cat << 'EOF' > src/main.rs
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead, Read};
use std::path::Path;
use flate2::read::GzDecoder;

#[derive(Debug, Deserialize, Serialize)]
struct Metadata {
    version: String,
    mutations: Vec<Mutation>,
    total_nodes: usize,
    config: Config,
}

#[derive(Debug, Deserialize, Serialize)]
#[serde(untagged)]
enum Mutation {
    AA {
        gene: String,
        previous_residue: String,
        residue_pos: usize,
        new_residue: String,
        mutation_id: usize,
        nuc_for_codon: usize,
        #[serde(rename = "type")]
        mutation_type: String,
    },
    NT {
        gene: String,
        previous_residue: String,
        residue_pos: usize,
        new_residue: String,
        mutation_id: usize,
        #[serde(rename = "type")]
        mutation_type: String,
    },
}

#[derive(Debug, Deserialize, Serialize)]
struct Config {
    gene_details: HashMap<String, GeneDetail>,
    num_tips: usize,
}

#[derive(Debug, Deserialize, Serialize)]
struct GeneDetail {
    name: String,
    strand: i32,
    start: usize,
    end: usize,
}

#[derive(Debug, Deserialize, Serialize)]
struct Node {
    name: String,
    x_dist: f64,
    y: f64,
    mutations: Vec<i32>,
    meta_genbank_accession: String,
    meta_date: String,
    meta_country: String,
    meta_pangolin_lineage: String,
    parent_id: i32,
    node_id: i32,
    num_tips: i32,
    clades: HashMap<String, String>,
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        println!("Usage: {} <path_to_jsonl_file>", args[0]);
        std::process::exit(1);
    }

    let path = Path::new(&args[1]);
    let file = File::open(&path)?;

    let mut reader: Box<dyn BufRead> = if path.extension().and_then(|s| s.to_str()) == Some("gz") {
        Box::new(io::BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(io::BufReader::new(file))
    };

    // Read the first line separately as metadata
    let mut metadata_line = String::new();
    reader.read_line(&mut metadata_line)?;
    let metadata: Metadata = serde_json::from_str(&metadata_line)?;

    println!("Metadata:");
    println!("Version: {}", metadata.version);
    println!("Total nodes: {}", metadata.total_nodes);
    println!("Number of tips: {}", metadata.config.num_tips);
    println!("Number of mutations: {}", metadata.mutations.len());

    // Count AA and NT mutations
    let (aa_mutations, nt_mutations): (Vec<&Mutation>, Vec<&Mutation>) = metadata.mutations.iter()
        .partition(|m| matches!(m, Mutation::AA { .. }));
    println!("Number of AA mutations: {}", aa_mutations.len());
    println!("Number of NT mutations: {}", nt_mutations.len());

    let mut nodes: Vec<Node> = Vec::new();

    // Read and parse the remaining lines
    for line in reader.lines() {
        let line = line?;
        let node: Node = serde_json::from_str(&line)?;
        nodes.push(node);
    }

    // Process the collected node data
    let total_nodes = nodes.len();
    let total_mutations: usize = nodes.iter().map(|node| node.mutations.len()).sum();
    let max_x_dist = nodes.iter().map(|node| node.x_dist).fold(f64::NEG_INFINITY, f64::max);
    let min_x_dist = nodes.iter().map(|node| node.x_dist).fold(f64::INFINITY, f64::min);

    // Find the root node (where parent_id == node_id)
    let root_node = nodes.iter().find(|node| node.parent_id == node.node_id);

    println!("\nNode Data Analysis:");
    println!("Total nodes: {}", total_nodes);
    println!("Total mutations in nodes: {}", total_mutations);
    println!("Max x_dist: {}", max_x_dist);
    println!("Min x_dist: {}", min_x_dist);

 

    // Example of additional processing using the collected nodes
    let avg_num_tips: f64 = nodes.iter().map(|node| node.num_tips as f64).sum::<f64>() / total_nodes as f64;
    println!("Average number of tips per node: {:.2}", avg_num_tips);

    Ok(())
}
EOF

# Build the project
cargo build --release

# Check if the tfci.jsonl or tfci.jsonl.gz file exists in the parent directory
if [ -f "../tfci.jsonl" ]; then
    echo "Running the JSONL processor on tfci.jsonl"
    ./target/release/jsonl_processor "../tfci.jsonl"
elif [ -f "../tfci.jsonl.gz" ]; then
    echo "Running the JSONL processor on tfci.jsonl.gz"
    ./target/release/jsonl_processor "../tfci.jsonl.gz"
else
    echo "Neither tfci.jsonl nor tfci.jsonl.gz found in the parent directory."
    echo "Usage: ./target/release/jsonl_processor <path_to_jsonl_file>"
fi

# Return to the original directory
cd ..

echo "Setup complete. You can now run the JSONL processor with:"
echo "cd $PROJECT_NAME && ./target/release/jsonl_processor <path_to_jsonl_file>"