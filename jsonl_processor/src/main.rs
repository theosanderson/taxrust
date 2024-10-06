use actix_web::{web, App, HttpServer, Responder, Result};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
use std::sync::Mutex;
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

#[derive(Debug, Deserialize, Serialize, Clone)]
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

struct AppState {
    nodes: Mutex<HashMap<i32, Node>>,
}

async fn get_node(data: web::Data<AppState>, node_id: web::Path<i32>) -> Result<impl Responder> {
    let nodes = data.nodes.lock().unwrap();
    if let Some(node) = nodes.get(&node_id) {
        Ok(web::Json(node.clone()))
    } else {
        Err(actix_web::error::ErrorNotFound("Node not found"))
    }
}

#[actix_web::main]
async fn main() -> std::io::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        println!("Usage: {} <path_to_jsonl_file>", args[0]);
        std::process::exit(1);
    }

    let path = Path::new(&args[1]);
    let file = File::open(&path).expect("Failed to open file");

    let reader: Box<dyn BufRead> = if path.extension().and_then(|s| s.to_str()) == Some("gz") {
        Box::new(io::BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(io::BufReader::new(file))
    };

    let mut lines = reader.lines();

    // Read the first line separately as metadata
    let metadata_line = lines.next().expect("Empty file").expect("Failed to read metadata");
    let _metadata: Metadata = serde_json::from_str(&metadata_line).expect("Failed to parse metadata");

    let mut nodes = HashMap::new();

    // Process nodes
    for line in lines {
        let line = line.expect("Failed to read line");
        let node: Node = serde_json::from_str(&line).expect("Failed to parse node");
        nodes.insert(node.node_id, node);
    }

    let app_state = web::Data::new(AppState {
        nodes: Mutex::new(nodes),
    });

    println!("Starting server at http://localhost:8080");

    HttpServer::new(move || {
        App::new()
            .app_data(app_state.clone())
            .route("/node/{node_id}", web::get().to(get_node))
    })
    .bind("127.0.0.1:8080")?
    .run()
    .await
}