## The sql command to update from EDGE v1.1 db schema to new db schema 

CREATE TABLE IF NOT EXISTS reports (
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

CREATE TABLE IF NOT EXISTS users_reports (
    user_id INT NOT NULL,
    report_id INT NOT NULL,
    role varchar(255) NOT NULL DEFAULT 'owner',
    PRIMARY KEY (user_id,report_id),
    FOREIGN KEY ur_user_id_key (user_id) REFERENCES users(id),
    FOREIGN KEY ur_report_id_key (report_id) REFERENCES reports(id)
);
ALTER TABLE users_reports ADD CONSTRAINT rroleChk CHECK (role IN ('owner','guest'));
ALTER TABLE reports ADD CONSTRAINT rpublishedChk CHECK (published IN ('yes','no'));

ALTER TABLE projects ADD full_name VARCHAR(255);
ALTER TABLE users_projects add display varchar(25) NOT NULL DEFAULT 'yes';
ALTER TABLE users_projects ADD CONSTRAINT displayChk CHECK (display IN ('yes','no'));

ALTER TABLE projects ADD run_submitted DATETIME;
ALTER TABLE projects ADD running_time TIME;

#clean up projects
DELETE FROM users_projects WHERE project_id in (select id FROM projects WHERE status='delete');
DELETE FROM projects WHERE status='delete';
