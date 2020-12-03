CREATE DATABASE IF NOT EXIST userManagement;
LOAD ~/edge/userManagement/userManagement_schema.sql;
LOAD ~/edge/userManagement/userManagement_constrains.sql;
CREATE USER 'edge'@'localhost' IDENTIFIED BY 'changePassword';
GRANT ALL PRIVILEGES ON userManagement.* TO 'edge'@'localhost';
