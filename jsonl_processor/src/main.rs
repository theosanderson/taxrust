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

#[derive(Debug, Deserialize)]
struct Node {
    name: String,
    x_dist: f64,
    y: f64,
    mutations: Vec<i32>,
    parent_id: i32,
    node_id: i32,
    num_tips: i32,
    clades: HashMap<String, String>,
    #[serde(flatten)]
    meta: HashMap<String, Value>,
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        println!("Usage: {} <path_to_jsonl_file>", args[0]);
        std::process::exit(1);
    }

    let path = Path::new(&args[1]);
    let file = File::open(&path)?;

    let reader: Box<dyn BufRead> = if path.extension().and_then(|s| s.to_str()) == Some("gz") {
        Box::new(io::BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(io::BufReader::new(file))
    };

    let mut lines = reader.lines();

    // Read the first line separately as metadata
    let metadata_line = lines.next().ok_or("Empty file")??;
    let metadata: Metadata = serde_json::from_str(&metadata_line)?;

    println!("Metadata:");
    println!("Version: {}", metadata.version);
    println!("Total nodes: {}", metadata.total_nodes);
    println!("Number of tips: {}", metadata.config.num_tips);
    println!("Number of mutations: {}", metadata.mutations.len());

    // Count AA and NT mutations
    let mut aa_mutations = 0;
    let mut nt_mutations = 0;
    for mutation in &metadata.mutations {
        match mutation {
            Mutation::AA { .. } => aa_mutations += 1,
            Mutation::NT { .. } => nt_mutations += 1,
        }
    }
    println!("Number of AA mutations: {}", aa_mutations);
    println!("Number of NT mutations: {}", nt_mutations);

    // Read the second line to determine dynamic fields
    let second_line = lines.next().ok_or("Missing second line")??;
    let sample_node: Node = serde_json::from_str(&second_line)?;
    let dynamic_fields: Vec<String> = sample_node.meta.keys()
        .filter(|&k| k.starts_with("meta_"))
        .cloned()
        .collect();

    println!("\nDynamic meta fields:");
    for field in &dynamic_fields {
        println!("  {}", field);
    }

    // Process nodes
    let mut total_nodes = 1; // Count the sample node
    let mut total_mutations = sample_node.mutations.len();
    let mut max_x_dist = sample_node.x_dist;
    let mut min_x_dist = sample_node.x_dist;
    let mut sum_num_tips = sample_node.num_tips as f64;
    let mut root_node_id = if sample_node.parent_id == sample_node.node_id {
        Some(sample_node.node_id)
    } else {
        None
    };

    // Process remaining lines
    for line in lines {
        let line = line?;
        let node: Node = serde_json::from_str(&line)?;
        
        total_nodes += 1;
        total_mutations += node.mutations.len();
        max_x_dist = max_x_dist.max(node.x_dist);
        min_x_dist = min_x_dist.min(node.x_dist);
        sum_num_tips += node.num_tips as f64;

        if node.parent_id == node.node_id {
            root_node_id = Some(node.node_id);
        }
    }

    println!("\nNode Data Analysis:");
    println!("Total nodes: {}", total_nodes);
    println!("Total mutations in nodes: {}", total_mutations);
    println!("Max x_dist: {}", max_x_dist);
    println!("Min x_dist: {}", min_x_dist);

    // Calculate average number of tips per node
    let avg_num_tips = sum_num_tips / total_nodes as f64;
    println!("Average number of tips per node: {:.2}", avg_num_tips);

    // Print sample data for dynamic fields
    if !dynamic_fields.is_empty() {
        println!("\nSample data for dynamic fields (from the second line):");
        for field in &dynamic_fields {
            if let Some(value) = sample_node.meta.get(field) {
                println!("  {}: {}", field, value);
            }
        }
    }

    Ok(())
}