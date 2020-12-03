DROP TABLE if EXISTS studies;
DROP TABLE if EXISTS study_types;
DROP TABLE if EXISTS animal_hosts;
DROP TABLE if EXISTS host_isolation_source;
DROP TABLE if EXISTS nonhost_isolation_source;
DROP TABLE if EXISTS seq_centers;
DROP TABLE if EXISTS sequencers;
DROP TABLE if EXISTS symptom_categories;
DROP TABLE if EXISTS symptoms;
DROP TABLE if EXISTS runs;


CREATE TABLE studies (
    id INT AUTO_INCREMENT PRIMARY KEY,
    name VARCHAR(765) NOT NULL UNIQUE
);
ALTER TABLE studies AUTO_INCREMENT=100001;

CREATE TABLE study_types (
    id INT AUTO_INCREMENT PRIMARY KEY,
    name VARCHAR(255) NOT NULL UNIQUE
);

CREATE TABLE animal_hosts (
    id INT AUTO_INCREMENT PRIMARY KEY,
    name VARCHAR(255) NOT NULL UNIQUE
);

CREATE TABLE host_isolation_source (
    id INT AUTO_INCREMENT PRIMARY KEY,
    name VARCHAR(255) NOT NULL UNIQUE
);

CREATE TABLE nonhost_isolation_source (
    id INT AUTO_INCREMENT PRIMARY KEY,
    name VARCHAR(255) NOT NULL UNIQUE
);

CREATE TABLE seq_centers (
    id INT AUTO_INCREMENT PRIMARY KEY,
    name VARCHAR(255) NOT NULL UNIQUE
);

CREATE TABLE sequencers (
    id INT AUTO_INCREMENT PRIMARY KEY,
    name VARCHAR(255) NOT NULL UNIQUE
);

CREATE TABLE symptom_categories (
    id INT NOT NULL PRIMARY KEY,
    name VARCHAR(255) NOT NULL UNIQUE
);

CREATE TABLE symptoms (
    id INT AUTO_INCREMENT PRIMARY KEY,
    cat_id INT NOT NULL,
    name VARCHAR(255) NOT NULL UNIQUE,
    FOREIGN KEY cat_id_key (cat_id) REFERENCES symptom_categories(id)
);

CREATE TABLE runs (
    id INT AUTO_INCREMENT PRIMARY KEY,
    name VARCHAR(255) NOT NULL
);
ALTER TABLE runs AUTO_INCREMENT=100001;

INSERT INTO study_types(name) VALUES('Whole Genome Sequencing');
INSERT INTO study_types(name) VALUES('Metagenomics');
INSERT INTO study_types(name) VALUES('Transcriptome Analysis');
INSERT INTO study_types(name) VALUES('Resequencing');
INSERT INTO study_types(name) VALUES('Epigenetics');
INSERT INTO study_types(name) VALUES('Synthetic Genomics');
INSERT INTO study_types(name) VALUES('Forensic or Paleo-genomics');
INSERT INTO study_types(name) VALUES('Gene Regulation Study');
INSERT INTO study_types(name) VALUES('Cancer Genomics');
INSERT INTO study_types(name) VALUES('Population Genomics');
INSERT INTO study_types(name) VALUES('RNASeq');
INSERT INTO study_types(name) VALUES('Exome Sequencing');
INSERT INTO study_types(name) VALUES('Pooled Clone Sequencing');

INSERT INTO animal_hosts(name) VALUES('Dog');

INSERT INTO host_isolation_source(name) VALUES('Blood');
INSERT INTO host_isolation_source(name) VALUES('Nasal');
INSERT INTO host_isolation_source(name) VALUES('Saliva (oral)');
INSERT INTO host_isolation_source(name) VALUES('Skin');
INSERT INTO host_isolation_source(name) VALUES('Sputum');
INSERT INTO host_isolation_source(name) VALUES('Stool (gut)');
INSERT INTO host_isolation_source(name) VALUES('Throat');
INSERT INTO host_isolation_source(name) VALUES('Vaginal');
INSERT INTO host_isolation_source(name) VALUES('Wound');
INSERT INTO host_isolation_source(name) VALUES('Unknown');

INSERT INTO nonhost_isolation_source(name) VALUES('Air');
INSERT INTO nonhost_isolation_source(name) VALUES('Built-environment');
INSERT INTO nonhost_isolation_source(name) VALUES('Microbial mat/biofilm');
INSERT INTO nonhost_isolation_source(name) VALUES('Plant');
INSERT INTO nonhost_isolation_source(name) VALUES('Sediment');
INSERT INTO nonhost_isolation_source(name) VALUES('Soil');
INSERT INTO nonhost_isolation_source(name) VALUES('Water');
INSERT INTO nonhost_isolation_source(name) VALUES('Wastewater/Sludge');
INSERT INTO nonhost_isolation_source(name) VALUES('Unknown');

INSERT INTO sequencers(name) VALUES('Illumina Genome Analyzer');
INSERT INTO sequencers(name) VALUES('Illumina HiSeq');
INSERT INTO sequencers(name) VALUES('Illumina MiSeq');
INSERT INTO sequencers(name) VALUES('Illumina NextSeq');
INSERT INTO sequencers(name) VALUES('Illumina MiniSeq');
INSERT INTO sequencers(name) VALUES('Ion Torrent S5');
INSERT INTO sequencers(name) VALUES('Ion Torrent PGM');
INSERT INTO sequencers(name) VALUES('Ion Torrent Proton');
INSERT INTO sequencers(name) VALUES('Nanopore MinION');
INSERT INTO sequencers(name) VALUES('PacBio RS');
INSERT INTO sequencers(name) VALUES('PacBio Sequel');

INSERT INTO symptom_categories(id,name) VALUES(1,'Constitutional');
INSERT INTO symptom_categories(id,name) VALUES(2,'Eyes');
INSERT INTO symptom_categories(id,name) VALUES(3,'ENT/Mouth');
INSERT INTO symptom_categories(id,name) VALUES(4,'Respiratory');
INSERT INTO symptom_categories(id,name) VALUES(5,'Cardiovascular');
INSERT INTO symptom_categories(id,name) VALUES(6,'Gastrointestinal');
INSERT INTO symptom_categories(id,name) VALUES(7,'Urinary');
INSERT INTO symptom_categories(id,name) VALUES(8,'Musculoskeletal');
INSERT INTO symptom_categories(id,name) VALUES(9,'Skin');
INSERT INTO symptom_categories(id,name) VALUES(10,'Neurologic');
INSERT INTO symptom_categories(id,name) VALUES(11,'Heme/Lymph');

INSERT INTO symptoms(cat_id,name) VALUES(1,'Recent weight loss or gain');
INSERT INTO symptoms(cat_id,name) VALUES(1,'Appetite changes');
INSERT INTO symptoms(cat_id,name) VALUES(1,'Fatigue');
INSERT INTO symptoms(cat_id,name) VALUES(1,'Fever');
INSERT INTO symptoms(cat_id,name) VALUES(1,'Shaking and chills');

INSERT INTO symptoms(cat_id,name) VALUES(2,'Eye redness');
INSERT INTO symptoms(cat_id,name) VALUES(2,'Eye drainage');
INSERT INTO symptoms(cat_id,name) VALUES(2,'Dry, irritated eyes');

INSERT INTO symptoms(cat_id,name) VALUES(3,'Ear pain');
INSERT INTO symptoms(cat_id,name) VALUES(3,'Sinus pain');
INSERT INTO symptoms(cat_id,name) VALUES(3,'Nosebleeds');
INSERT INTO symptoms(cat_id,name) VALUES(3,'Dizziness');
INSERT INTO symptoms(cat_id,name) VALUES(3,'Sore throat');
INSERT INTO symptoms(cat_id,name) VALUES(3,'Hoarseness');
INSERT INTO symptoms(cat_id,name) VALUES(3,'Difficulty swallowing');

INSERT INTO symptoms(cat_id,name) VALUES(4,'Blood in sputum');
INSERT INTO symptoms(cat_id,name) VALUES(4,'Chest tightness or pain');
INSERT INTO symptoms(cat_id,name) VALUES(4,'Cough');
INSERT INTO symptoms(cat_id,name) VALUES(4,'Shortness of breath');
INSERT INTO symptoms(cat_id,name) VALUES(4,'Wheezing');

INSERT INTO symptoms(cat_id,name) VALUES(5,'Chest pain or heaviness');
INSERT INTO symptoms(cat_id,name) VALUES(5,'Swelling of feet or legs');

INSERT INTO symptoms(cat_id,name) VALUES(6,'Abdominal pain or tenderness');
INSERT INTO symptoms(cat_id,name) VALUES(6,'Blood in stool');
INSERT INTO symptoms(cat_id,name) VALUES(6,'Diarrhea');
INSERT INTO symptoms(cat_id,name) VALUES(6,'Vomiting');
INSERT INTO symptoms(cat_id,name) VALUES(6,'Nausea');

INSERT INTO symptoms(cat_id,name) VALUES(7,'Blood in urine');
INSERT INTO symptoms(cat_id,name) VALUES(7,'Painful or difficult urination');

INSERT INTO symptoms(cat_id,name) VALUES(8,'Joint pain or swelling');
INSERT INTO symptoms(cat_id,name) VALUES(8,'Muscle aches');
INSERT INTO symptoms(cat_id,name) VALUES(8,'Muscle weakness');
INSERT INTO symptoms(cat_id,name) VALUES(8,'Neck pain');
INSERT INTO symptoms(cat_id,name) VALUES(8,'Back pain');

INSERT INTO symptoms(cat_id,name) VALUES(9,'Rash');

INSERT INTO symptoms(cat_id,name) VALUES(10,'Sezures');
INSERT INTO symptoms(cat_id,name) VALUES(10,'Headache');
INSERT INTO symptoms(cat_id,name) VALUES(10,'Extremity pain, numbness, tingling');
INSERT INTO symptoms(cat_id,name) VALUES(10,'Difficulty falling or staying asleep');

INSERT INTO symptoms(cat_id,name) VALUES(11,'Gum bleeding');
INSERT INTO symptoms(cat_id,name) VALUES(11,'Unexplained bruishing');
INSERT INTO symptoms(cat_id,name) VALUES(11,'Swollen, painful lymph nodes');


