DROP TABLE IF EXISTS node;
DROP TABLE IF EXISTS edge;
DROP TABLE IF EXISTS node_identifier;
DROP TABLE IF EXISTS tissue;
DROP TABLE IF EXISTS node_tissue;


CREATE TABLE node (
    id INTEGER PRIMARY KEY,
    name TEXT NOT NULL,
    display_name TEXT NOT NULL,
    tax_id INTEGER,
    type TEXT NOT NULL DEFAULT 'protein',
    pathways TEXT,
    source TEXT,
    function TEXT,
    source_database TEXT
);

CREATE TABLE tissue (
    id INTEGER PRIMARY KEY,
    tissue_id TEXT NOT NULL UNIQUE,
    ontology TEXT NOT NULL,
    name TEXT NOT NULL,
    full_name TEXT
);

CREATE TABLE node_tissue (
    node_id INTEGER NOT NULL,
    tissue_id INTEGER NOT NULL,
    PRIMARY KEY (node_id, tissue_id),
    FOREIGN KEY (node_id) REFERENCES node(id),
    FOREIGN KEY (tissue_id) REFERENCES tissue(id)
);

CREATE TABLE edge (
    id INTEGER PRIMARY KEY,
    interactor_a_node_id INTEGER NOT NULL,
    interactor_b_node_id INTEGER NOT NULL,
    interactor_a_node_name TEXT NOT NULL,
    interactor_b_node_name TEXT NOT NULL,
    layer INTEGER NOT NULL,
    interaction_types TEXT,
    source_db TEXT NOT NULL,
    pathway_type TEXT,
    is_directed BOOLEAN DEFAULT 0,
    publications TEXT,
    FOREIGN KEY(interactor_a_node_id) REFERENCES node (id),
    FOREIGN KEY(interactor_b_node_id) REFERENCES node (id)
);

CREATE INDEX idx_node_name ON node(name);
CREATE INDEX idx_node_name_tax ON node(name, tax_id);
CREATE INDEX idx_node_type ON node(type);
CREATE INDEX idx_node_source ON node(source_database);
CREATE INDEX idx_node_type_source ON node(type, source_database);
CREATE INDEX idx_tissue_ontology ON tissue(ontology);
CREATE INDEX idx_tissue_id ON tissue(tissue_id);
CREATE INDEX idx_node_tissue_tissue ON node_tissue(tissue_id);
CREATE INDEX idx_node_tissue_node ON node_tissue(node_id);
CREATE INDEX idx_edge_interactors ON edge(interactor_a_node_id, interactor_b_node_id);
CREATE INDEX idx_edge_layer ON edge(layer);
CREATE INDEX idx_edge_pathway ON edge(pathway_type);
CREATE INDEX idx_edge_source_pathway ON edge(source_db, pathway_type);
