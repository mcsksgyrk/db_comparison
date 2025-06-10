DROP TABLE IF EXISTS `node`;
DROP TABLE IF EXISTS `edge`;
DROP TABLE IF EXISTS `node_identifier`;
DROP TABLE IF EXISTS `tissue`;
DROP TABLE IF EXISTS `node_tissue`;

PRAGMA foreign_keys = ON;

-- Core node table
CREATE TABLE `node` (
    `id` INTEGER PRIMARY KEY,
    `name` TEXT NOT NULL,         -- Primary identifier (UniProt/PubChem/KEGG ID)
    `primary_id_type` TEXT,       -- Type of the primary identifier in name field
    `tax_id` INTEGER,             -- Only for biological entities, now nullable
    `type` TEXT NOT NULL DEFAULT 'protein', -- 'protein', 'small_molecule', 'rna', etc.
    `pathways` TEXT,              -- JSON array of pathway IDs
    `source` TEXT,                -- Original data source
    `function` TEXT,              -- Functional description
    `source_database` TEXT        -- 'autophagy', 'ferroptosis', 'both'
);

CREATE TABLE `node_identifier` (
    `node_id` INTEGER NOT NULL,
    `id_type` TEXT NOT NULL,      -- 'kegg_id', 'uniprot_id', 'pubchem_id', etc.
    `is_primary` BOOLEAN DEFAULT 0, -- Indicates if this is the primary ID used in node.name
    `id_value` TEXT NOT NULL,
    PRIMARY KEY (`node_id`, `id_type`),
    FOREIGN KEY (`node_id`) REFERENCES `node`(`id`) ON UPDATE NO ACTION ON DELETE CASCADE
);

CREATE TABLE `tissue` (
    `id` INTEGER PRIMARY KEY,
    `tissue_id` TEXT NOT NULL UNIQUE,     -- 'UBERON:0000955', 'CL:0000738'
    `ontology` TEXT NOT NULL,             -- 'UBERON', 'CL'
    `name` TEXT NOT NULL,                 -- 'brain', 'leukocyte'
    `full_name` TEXT                      -- Cleaned full name
);

CREATE TABLE `node_tissue` (
    `node_id` INTEGER NOT NULL,
    `tissue_id` INTEGER NOT NULL,
    PRIMARY KEY (`node_id`, `tissue_id`),
    FOREIGN KEY (`node_id`) REFERENCES `node`(`id`) ON DELETE CASCADE,
    FOREIGN KEY (`tissue_id`) REFERENCES `tissue`(`id`) ON DELETE CASCADE
);

-- Interaction edges
CREATE TABLE `edge` (
    `id` INTEGER PRIMARY KEY,
    `interactor_a_node_id` INTEGER NOT NULL,
    `interactor_b_node_id` INTEGER NOT NULL,
    `interactor_a_node_name` TEXT NOT NULL,
    `interactor_b_node_name` TEXT NOT NULL,
    `layer` TEXT NOT NULL,
    `interaction_types` TEXT,     -- JSON array of interaction types
    `source_db` TEXT NOT NULL,
    `pathway_type` TEXT,          -- 'autophagy', 'ferroptosis', 'cross-pathway'
    `confidence_score` REAL,      -- Interaction confidence score
    `is_directed` BOOLEAN DEFAULT 0, -- Whether interaction is directional
    `publications` TEXT,          -- JSON array of publication references
    FOREIGN KEY(`interactor_a_node_id`) REFERENCES node (`id`) ON UPDATE NO ACTION ON DELETE CASCADE,
    FOREIGN KEY(`interactor_b_node_id`) REFERENCES node (`id`) ON UPDATE NO ACTION ON DELETE CASCADE
);

-- Indexes for performance
CREATE INDEX `idx_node_name` ON `node`(`name`);
CREATE INDEX `idx_node_name_tax` ON `node`(`name`, `tax_id`);
CREATE INDEX `idx_node_type` ON `node`(`type`);
CREATE INDEX `idx_node_source` ON `node`(`source_database`);
CREATE INDEX `idx_node_type_source` ON `node`(`type`, `source_database`);

CREATE INDEX `idx_node_identifier_type` ON `node_identifier`(`id_type`);
CREATE INDEX `idx_node_identifier_value` ON `node_identifier`(`id_value`);
CREATE INDEX `idx_node_identifier_node` ON `node_identifier`(`node_id`);

CREATE INDEX `idx_tissue_ontology` ON `tissue`(`ontology`);
CREATE INDEX `idx_tissue_id` ON `tissue`(`tissue_id`);

CREATE INDEX `idx_node_tissue_tissue` ON `node_tissue`(`tissue_id`);
CREATE INDEX `idx_node_tissue_node` ON `node_tissue`(`node_id`);

CREATE INDEX `idx_edge_interactors` ON `edge`(`interactor_a_node_id`, `interactor_b_node_id`);
CREATE INDEX `idx_edge_layer` ON `edge`(`layer`);
CREATE INDEX `idx_edge_pathway` ON `edge`(`pathway_type`);
CREATE INDEX `idx_edge_source_pathway` ON `edge`(`source_db`, `pathway_type`);
