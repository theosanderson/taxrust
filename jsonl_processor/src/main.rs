use actix_web::{web, App, HttpServer, Responder, Result, get, HttpResponse};
use actix_cors::Cors;
use serde::{Deserialize, Serialize};
use serde_json::json;
use serde_json::Value;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::io::Seek;
use std::io::{self, BufRead};
use std::path::Path;
use std::sync::Mutex;
use flate2::read::GzDecoder;
use std::cmp::Ordering;
use std::time::Instant;

#[derive(Debug, Deserialize, Serialize, Clone)]
struct Metadata {
    version: String,
    mutations: Vec<Mutation>,
    total_nodes: usize,
    config: Config,
}

#[derive(Debug, Deserialize, Serialize, Clone)]
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

#[derive(Debug, Deserialize, Serialize, Clone)]
struct Config {
    gene_details: HashMap<String, GeneDetail>,
    num_tips: usize,
    #[serde(default)]
    mutations: Vec<Mutation>,
    #[serde(default)]
    initial_x: Option<f64>,
    #[serde(default)]
    initial_y: Option<f64>,
    #[serde(default)]
    initial_zoom: Option<f64>,
    #[serde(default)]
    keys_to_display: Option<Vec<String>>,
    #[serde(default)]
    num_nodes: Option<usize>,
    #[serde(default)]
    root_mutations: Option<Vec<i32>>,
    #[serde(default)]
    root_id: Option<i32>,
}

#[derive(Debug, Deserialize, Serialize, Clone)]
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
    child_to_parent: Mutex<HashMap<i32, i32>>,
    config: Mutex<Config>,
    root_mutations: Mutex<Vec<i32>>,
    root_id: Mutex<i32>,
}

#[derive(Debug, Deserialize)]
struct NodesQuery {
    min_y: Option<f64>,
    max_y: Option<f64>,
    min_x: Option<f64>,
    max_x: Option<f64>,
    x_type: Option<String>,
}

#[derive(Debug, Serialize)]
struct NodesResponse {
    nodes: Vec<Node>,
}

fn load_data(path: &Path) -> Result<(Metadata, HashMap<i32, Node>, HashMap<i32, i32>, Vec<i32>, i32), Box<dyn Error>> {
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
    let mut child_to_parent = HashMap::new();
    let mut root_mutations = Vec::new();
    let mut root_id = 0;

    // Process nodes
    for line in lines {
        let line = line?;
        let mut node: Node = serde_json::from_str(&line)?;
        
        if node.parent_id == node.node_id {
            // This is the root node
            root_mutations = node.mutations.clone();
            root_id = node.node_id;
            node.mutations = Vec::new(); // Clear root node mutations
        } else {
            child_to_parent.insert(node.node_id, node.parent_id);
        }
        
        nodes.insert(node.node_id, node);
    }

    Ok((metadata, nodes, child_to_parent, root_mutations, root_id))
}

fn scale_y_coordinates(nodes: &mut HashMap<i32, Node>) {
    let num_nodes = nodes.len();
    let scale_y = 24e2 / if num_nodes > 10000 { num_nodes as f64 } else { num_nodes as f64 * 0.6666 };
    
    for node in nodes.values_mut() {
        node.y = (node.y * scale_y * 1e6).round() / 1e6;  // Round to 6 decimal places
    }
}

fn calculate_extremes(nodes: &HashMap<i32, Node>) -> (f64, f64, f64, f64) {
    let mut min_y = f64::MAX;
    let mut max_y = f64::MIN;
    let mut min_x = f64::MAX;
    let mut max_x = f64::MIN;

    for node in nodes.values() {
        min_y = min_y.min(node.y);
        max_y = max_y.max(node.y);
        min_x = min_x.min(node.x_dist);
        max_x = max_x.max(node.x_dist);
    }

    (min_y, max_y, min_x, max_x)
}

fn update_config(config: &mut Config, nodes: &HashMap<i32, Node>, root_mutations: &Vec<i32>, root_id: i32, mutations: Vec<Mutation>) {
    let (min_y, max_y, min_x, max_x) = calculate_extremes(nodes);
    config.initial_x = Some((max_x + min_x) / 2.0);
    config.initial_y = Some((max_y + min_y) / 2.0);
    config.initial_zoom = Some(config.initial_zoom.unwrap_or(-2.0));
    config.num_nodes = Some(nodes.len());
    config.root_mutations = Some(root_mutations.clone());
    config.root_id = Some(root_id);
    config.mutations = mutations;
    config.keys_to_display = Some(vec!["name".to_string(), "num_tips".to_string()]);
}

#[get("/config/")]
async fn get_config(data: web::Data<AppState>) -> impl Responder {
    let config = data.config.lock().unwrap();
    HttpResponse::Ok().json(config.clone())
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

#[get("/search/")]
async fn search(_data: web::Data<AppState>) -> impl Responder {
    HttpResponse::Ok().json(json!({
        "type": "complete",
        "data": [],
        "total_count": 0
    }))
}

#[get("/nodes/")]
async fn get_nodes(
    data: web::Data<AppState>,
    query: web::Query<NodesQuery>,
) -> impl Responder {
    let start_time = Instant::now();

    let nodes = data.nodes.lock().unwrap();
    let child_to_parent = data.child_to_parent.lock().unwrap();
    let lock_time = start_time.elapsed();
    println!("Time to acquire locks: {:?}", lock_time);
    
    let config = data.config.lock().unwrap();
    
    let min_y = query.min_y.unwrap_or_else(|| nodes.values().map(|n| n.y).min_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal)).unwrap_or(0.0));
    let max_y = query.max_y.unwrap_or_else(|| nodes.values().map(|n| n.y).max_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal)).unwrap_or(0.0));
    let min_x = query.min_x.unwrap_or_else(|| nodes.values().map(|n| n.x_dist).min_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal)).unwrap_or(0.0));
    let max_x = query.max_x.unwrap_or_else(|| nodes.values().map(|n| n.x_dist).max_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal)).unwrap_or(0.0));
    let x_type = query.x_type.as_deref().unwrap_or("x_dist");

    let query_time = start_time.elapsed() - lock_time;
    println!("Time to process query parameters: {:?}", query_time);

    println!("min_y: {}, max_y: {}, min_x: {}, max_x: {}", min_y, max_y, min_x, max_x);

    let filter_start = Instant::now();
    let filtered = filter_nodes(&nodes, min_y, max_y);
    let filter_time = filter_start.elapsed();
    println!("Time to filter nodes: {:?}", filter_time);

    let reduce_start = Instant::now();
    let reduced_leaves = reduce_overplotting(
        filtered.into_iter().filter(|n| n.num_tips == 1).collect(),
        get_precision(min_x, max_x),
        get_precision(min_y, max_y),
        x_type,
    );
    let reduce_time = reduce_start.elapsed();
    println!("Time to reduce overplotting: {:?}", reduce_time);

    let parents_start = Instant::now();
    let result = add_parents(&nodes, &child_to_parent, reduced_leaves);
    let parents_time = parents_start.elapsed();
    println!("Time to add parents: {:?}", parents_time);

    let total_time = start_time.elapsed();
    println!("Total time for /nodes/ endpoint: {:?}", total_time);

    HttpResponse::Ok().json(NodesResponse { nodes: result })
}

fn filter_nodes<'a>(nodes: &'a HashMap<i32, Node>, min_y: f64, max_y: f64) -> Vec<&'a Node> {
    nodes.values()
        .filter(|n| n.y >= min_y && n.y <= max_y)
        .collect()
}

fn get_precision(min: f64, max: f64) -> f64 {
    2000.0 / (max - min)
}

fn reduce_overplotting<'a>(nodes: Vec<&'a Node>, precision_x: f64, precision_y: f64, x_type: &str) -> Vec<&'a Node> {
    println!("Precision: {}, {}", precision_x, precision_y);
    println!("Before: {}", nodes.len());
    let precision_x = precision_x / 5.0;
    let mut included_points = HashMap::new();
    let result: Vec<&Node> = nodes.into_iter().filter(|node| {
        let x = if x_type == "x_dist" { node.x_dist } else { node.x_dist };
        let rounded_x = (x * precision_x).round() as i64;
        let rounded_y = (node.y * precision_y).round() as i64;
        included_points
            .entry(rounded_x)
            .or_insert_with(HashSet::new)
            .insert(rounded_y)
    }).collect();
    println!("After: {}", result.len());
    result
}

fn add_parents(all_nodes: &HashMap<i32, Node>, child_to_parent: &HashMap<i32, i32>, filtered: Vec<&Node>) -> Vec<Node> {
    let start = Instant::now();
    
    println!("Before: {}", filtered.len());
    let mut selected_node_ids: HashSet<i32> = filtered.iter().map(|n| n.node_id).collect();
    let starting_size = selected_node_ids.len();

    let setup_time = start.elapsed();
    println!("Setup time: {:?}", setup_time);

    let processing_start = Instant::now();
    let mut to_process: Vec<i32> = selected_node_ids.iter().cloned().collect();

    while let Some(node_id) = to_process.pop() {
        if let Some(&parent_id) = child_to_parent.get(&node_id) {
            if !selected_node_ids.contains(&parent_id) {
                selected_node_ids.insert(parent_id);
                to_process.push(parent_id);
            }
        }
    }
    let processing_time = processing_start.elapsed();
    println!("Parent processing time: {:?}", processing_time);

    let result_start = Instant::now();
    let result: Vec<Node> = selected_node_ids.iter()
        .filter_map(|node_id| all_nodes.get(node_id).cloned())
        .collect();
    let result_time = result_start.elapsed();
    println!("Result collection time: {:?}", result_time);

    println!("After: {}", result.len());
    println!("Went from {} to {} nodes.", starting_size, result.len());
    
    let total_time = start.elapsed();
    println!("Total execution time: {:?}", total_time);
    
    result
}

#[actix_web::main]
async fn main() -> std::io::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        println!("Usage: {} <path_to_jsonl_file>", args[0]);
        std::process::exit(1);
    }

    let path = Path::new(&args[1]);
    
    let (mut metadata, mut nodes, child_to_parent, root_mutations, root_id) = load_data(path).expect("Failed to load data");

    scale_y_coordinates(&mut nodes);
    update_config(&mut metadata.config, &nodes, &root_mutations, root_id, metadata.mutations.clone());
    let app_state = web::Data::new(AppState {
        nodes: Mutex::new(nodes),
        child_to_parent: Mutex::new(child_to_parent),
        config: Mutex::new(metadata.config),
        root_mutations: Mutex::new(root_mutations),
        root_id: Mutex::new(root_id),
    });

    println!("Starting server at http://localhost:8080");

    HttpServer::new(move || {
        // Create a CORS middleware
        let cors = Cors::default()
            .allow_any_origin()
            .allow_any_method()
            .allow_any_header()
            .max_age(3600);

        App::new()
            .wrap(cors)  // Wrap the entire application with CORS middleware
            .app_data(app_state.clone())
            .service(index)
            .service(get_node)
            .service(get_nodes)
            .service(get_config)
            .service(search)
    })
    .bind(("127.0.0.1", 8080))?
    .disable_signals()
    .run()
    .await
}