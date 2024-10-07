use actix_web::{web, App, HttpServer, Responder, Result, get};
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

fn load_data(path: &Path) -> Result<(Metadata, HashMap<i32, Node>), Box<dyn Error>> {
    let file = File::open(path)?;

    let reader: Box<dyn BufRead> = if path.extension().and_then(|s| s.to_str()) == Some("gz") {
        Box::new(io::BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(io::BufReader::new(file))
    };

    let mut lines = reader.lines();

    // Read the first line separately as metadata
    let metadata_line = lines.next().ok_or("Empty file")??;
    let metadata: Metadata = serde_json::from_str(&metadata_line)?;

    let mut nodes = HashMap::new();

    // Process nodes
    for line in lines {
        let line = line?;
        let node: Node = serde_json::from_str(&line)?;
        nodes.insert(node.node_id, node);
    }

    Ok((metadata, nodes))
}

#[get("/node/{node_id}")]
async fn get_node(data: web::Data<AppState>, node_id: web::Path<i32>) -> Result<impl Responder> {
    let nodes = data.nodes.lock().unwrap();
    if let Some(node) = nodes.get(&node_id) {
        Ok(web::Json(node.clone()))
    } else {
        Err(actix_web::error::ErrorNotFound("Node not found"))
    }
}

#[get("/")]
async fn index(_data: web::Data<AppState>) -> String {
    "Hello world!".to_string()
}

#[actix_web::main]
async fn main() -> std::io::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        println!("Usage: {} <path_to_jsonl_file>", args[0]);
        std::process::exit(1);
    }

    let path = Path::new(&args[1]);
    
    let (_metadata, nodes) = load_data(path).expect("Failed to load data");

    let app_state = web::Data::new(AppState {
        nodes: Mutex::new(nodes),
    });

    println!("Starting server at http://localhost:8080");

    HttpServer::new(move || {
        App::new()
            .app_data(app_state.clone())
            .service(index)
            .service(get_node)
    })
    .bind(("127.0.0.1", 8080))?
    .disable_signals()
    .run()
    .await
}