DROP TABLE if EXISTS users_projects;
DROP TABLE if EXISTS users_reports;
DROP TABLE if EXISTS users;
DROP TABLE if EXISTS projects;
DROP TABLE if EXISTS reports;

CREATE TABLE users (
    id INT AUTO_INCREMENT PRIMARY KEY,
    email VARCHAR(255) NOT NULL UNIQUE,
    password VARCHAR(255) NOT NULL,
    spassword VARCHAR(255) NOT NULL,
    first_name VARCHAR(255),
    last_name VARCHAR(255),
    type varchar(255) NOT NULL DEFAULT 'user',
    active varchar(25) NOT NULL DEFAULT 'no',
    created DATETIME,
    updated DATETIME
);

CREATE TABLE projects (
    id INT AUTO_INCREMENT PRIMARY KEY,
    name VARCHAR(255) NOT NULL,
    description VARCHAR(3000),
    status varchar(255) NOT NULL DEFAULT 'in process',
    published varchar(25) NOT NULL DEFAULT 'no',
    code VARCHAR(255) NOT NULL,
	full_name VARCHAR(255),
	run_submitted DATETIME,
	running_time TIME,
    created DATETIME,
    updated DATETIME
);

CREATE TABLE users_projects (
    user_id INT NOT NULL,
    project_id INT NOT NULL,
    role varchar(255) NOT NULL DEFAULT 'owner',
    display varchar(25) NOT NULL DEFAULT 'yes',
    PRIMARY KEY (user_id,project_id),
    FOREIGN KEY up_user_id_key (user_id) REFERENCES users(id),
    FOREIGN KEY up_project_id_key (project_id) REFERENCES projects(id)
);

CREATE TABLE reports (
    id INT AUTO_INCREMENT PRIMARY KEY,
    name VARCHAR(255) NOT NULL,
    description VARCHAR(3000),
    status varchar(255) NOT NULL DEFAULT 'in process',
    published varchar(25) NOT NULL DEFAULT 'no',
    code VARCHAR(255) NOT NULL,
    full_name VARCHAR(255),
    created DATETIME,
    updated DATETIME
);

CREATE TABLE users_reports (
    user_id INT NOT NULL,
    report_id INT NOT NULL,
    role varchar(255) NOT NULL DEFAULT 'owner',
    PRIMARY KEY (user_id,report_id),
    FOREIGN KEY ur_user_id_key (user_id) REFERENCES users(id),
    FOREIGN KEY ur_report_id_key (report_id) REFERENCES reports(id)
);

