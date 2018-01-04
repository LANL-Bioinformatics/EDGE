##constrains
ALTER TABLE users ADD CONSTRAINT utypeChk CHECK (type IN ('admin','user'));
ALTER TABLE users ADD CONSTRAINT actChk CHECK (active IN ('yes','no'));
ALTER TABLE users_projects ADD CONSTRAINT roleChk CHECK (role IN ('owner','guest'));
ALTER TABLE users_projects ADD CONSTRAINT displayChk CHECK (display IN ('yes','no'));
ALTER TABLE users_reports ADD CONSTRAINT rroleChk CHECK (role IN ('owner','guest'));
ALTER TABLE projects AUTO_INCREMENT = 100000;
ALTER TABLE projects ADD CONSTRAINT publishedChk CHECK (published IN ('yes','no'));
ALTER TABLE reports ADD CONSTRAINT rpublishedChk CHECK (published IN ('yes','no'));

