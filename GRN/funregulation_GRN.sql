--
-- File generated with SQLiteStudio v3.2.1 on seg mai 4 08:47:27 2020
--
-- Text encoding used: UTF-8
--
PRAGMA foreign_keys = off;
BEGIN TRANSACTION;

-- Table: gene
CREATE TABLE gene (
    organism  TEXT    NOT NULL,
    locus_tag TEXT    NOT NULL,
    name      TEXT    NOT NULL,
    is_tf     INTEGER DEFAULT 0,
    PRIMARY KEY (
        organism,
        locus_tag
    )
);


-- Table: network_node
CREATE TABLE network_node (
    organism       TEXT      NOT NULL,
    tf_locus_tag   TEXT      NOT NULL,
    tf_name        TEXT      NOT NULL,
    tf_ortho_names TEXT (20),
    tg_locus_tag   TEXT      NOT NULL,
    tg_name        TEXT      NOT NULL,
    tg_ortho_names TEXT (20),
    interaction    TEXT,
    ortho_pox      INTEGER,
    ortho_ani      INTEGER,
    ortho_ncr      INTEGER,
    ortho_sce      INTEGER,
    tfbs_count     INTEGER,
    PRIMARY KEY (
        organism,
        tf_locus_tag,
        tg_locus_tag
    ),
    FOREIGN KEY (
        organism
    )
    REFERENCES regulation (src_organism) ON DELETE RESTRICT
                                         ON UPDATE NO ACTION,
    FOREIGN KEY (
        tf_locus_tag
    )
    REFERENCES regulation (src_tf_locus_tag) ON DELETE RESTRICT
                                             ON UPDATE NO ACTION,
    FOREIGN KEY (
        tg_locus_tag
    )
    REFERENCES regulation (src_tg_locus_tag) ON DELETE RESTRICT
                                             ON UPDATE NO ACTION
);


-- Table: ortho
CREATE TABLE ortho (
    src_organism    TEXT NOT NULL,
    src_locus_tag   TEXT NOT NULL,
    ortho_organism  TEXT NOT NULL,
    ortho_locus_tag TEXT NOT NULL,
    PRIMARY KEY (
        src_organism,
        src_locus_tag,
        ortho_organism,
        ortho_locus_tag
    ),
    FOREIGN KEY (
        src_organism
    )
    REFERENCES gene (organism) ON DELETE RESTRICT
                               ON UPDATE NO ACTION,
    FOREIGN KEY (
        src_locus_tag
    )
    REFERENCES gene (locus_tag) ON DELETE RESTRICT
                                ON UPDATE NO ACTION,
    FOREIGN KEY (
        ortho_organism
    )
    REFERENCES gene (organism) ON DELETE RESTRICT
                               ON UPDATE NO ACTION,
    FOREIGN KEY (
        ortho_locus_tag
    )
    REFERENCES gene (locus_tag) ON DELETE RESTRICT
                                ON UPDATE NO ACTION
);


-- Table: pwm
CREATE TABLE pwm (
    organism       TEXT NOT NULL,
    locus_tag      TEXT NOT NULL,
    motif_id       TEXT NOT NULL,
    tf_status      TEXT,
    family_name    TEXT,
    motif_type     TEXT,
    msource_author TEXT,
    msource_year   TEXT,
    pmid           TEXT,
    PRIMARY KEY (
        organism,
        locus_tag,
        motif_id
    ),
    FOREIGN KEY (
        organism
    )
    REFERENCES gene (organism) ON DELETE RESTRICT
                               ON UPDATE NO ACTION,
    FOREIGN KEY (
        locus_tag
    )
    REFERENCES gene (locus_tag) ON DELETE RESTRICT
                                ON UPDATE NO ACTION
);


-- Table: regulation
CREATE TABLE regulation (
    src_organism       TEXT NOT NULL,
    src_tf_locus_tag   TEXT NOT NULL,
    src_tg_locus_tag   TEXT NOT NULL,
    interaction        TEXT,
    ortho_organism     TEXT NOT NULL,
    ortho_tf_locus_tag TEXT NOT NULL,
    ortho_tg_locus_tag TEXT NOT NULL,
    PRIMARY KEY (
        src_organism,
        src_tf_locus_tag,
        src_tg_locus_tag,
        ortho_organism,
        ortho_tf_locus_tag,
        ortho_tg_locus_tag
    ),
    FOREIGN KEY (
        src_organism
    )
    REFERENCES ortho (src_organism) ON DELETE RESTRICT
                                    ON UPDATE NO ACTION,
    FOREIGN KEY (
        src_tf_locus_tag
    )
    REFERENCES ortho (src_tf_locus_tag) ON DELETE RESTRICT
                                        ON UPDATE NO ACTION,
    FOREIGN KEY (
        src_tg_locus_tag
    )
    REFERENCES ortho (src_tg_locus_tag) ON DELETE RESTRICT
                                        ON UPDATE NO ACTION,
    FOREIGN KEY (
        ortho_organism
    )
    REFERENCES ortho (src_organism) ON DELETE RESTRICT
                                    ON UPDATE NO ACTION,
    FOREIGN KEY (
        ortho_tf_locus_tag
    )
    REFERENCES ortho (src_tf_locus_tag) ON DELETE RESTRICT
                                        ON UPDATE NO ACTION,
    FOREIGN KEY (
        ortho_tg_locus_tag
    )
    REFERENCES ortho (src_tg_locus_tag) ON DELETE RESTRICT
                                        ON UPDATE NO ACTION
);


-- Table: tfbs_prediction
CREATE TABLE tfbs_prediction (
    tfbs_id            INTEGER UNIQUE
                               PRIMARY KEY ASC AUTOINCREMENT,
    src_organism       TEXT    NOT NULL,
    src_tf_locus_tag   TEXT    NOT NULL,
    src_tg_locus_tag   TEXT    NOT NULL,
    ortho_organism     TEXT    NOT NULL,
    ortho_tf_locus_tag TEXT    NOT NULL,
    motif_id           TEXT,
    strand             TEXT,
    start              INTEGER,
    [end]              INTEGER,
    sequence           TEXT,
    weight             NUMBER,
    pval               TEXT,
    ln_pval            NUMBER,
    sig                NUMBER,
    FOREIGN KEY (
        src_organism
    )
    REFERENCES regulation (src_organism) ON DELETE NO ACTION
                                         ON UPDATE NO ACTION,
    FOREIGN KEY (
        src_organism
    )
    REFERENCES pwm (organism) ON DELETE NO ACTION
                              ON UPDATE NO ACTION,
    FOREIGN KEY (
        src_tf_locus_tag
    )
    REFERENCES regulation (src_tf_locus_tag) ON DELETE NO ACTION
                                             ON UPDATE NO ACTION,
    FOREIGN KEY (
        src_tg_locus_tag
    )
    REFERENCES regulation (src_tg_locus_tag) ON DELETE NO ACTION
                                             ON UPDATE NO ACTION,
    FOREIGN KEY (
        src_tf_locus_tag
    )
    REFERENCES pwm (locus_tag) ON DELETE NO ACTION
                               ON UPDATE NO ACTION,
    FOREIGN KEY (
        ortho_organism
    )
    REFERENCES regulation (ortho_organism) ON DELETE NO ACTION
                                           ON UPDATE NO ACTION,
    FOREIGN KEY (
        ortho_tf_locus_tag
    )
    REFERENCES regulation (ortho_tf_locus_tag) ON DELETE NO ACTION
                                               ON UPDATE NO ACTION,
    FOREIGN KEY (
        motif_id
    )
    REFERENCES pwm (motif_id) ON DELETE NO ACTION
                              ON UPDATE NO ACTION
);


COMMIT TRANSACTION;
PRAGMA foreign_keys = on;
