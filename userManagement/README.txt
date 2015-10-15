#install sendmail
sudo apt-get install sendmail

#setup mysql database
sudo apt-get install mysql-client mysql-server phpmyadmin

1. create database: userManagement
%mysql -p -u root
mysql> create database userManagement;
mysql> use userManagement;

2. load userManagement_schema.sql
    
3. load userManagement_constrains.sql

4. create an user account: 
              username: userManagement
              password: MtMmG3Ds2dzyYeQ5
   and grant all privileges on database "userManagement" to user userManagement
   
mysql> CREATE USER 'userManagement'@'localhost' IDENTIFIED BY 'MtMmG3Ds2dzyYeQ5';
Query OK, 0 rows affected (0.00 sec)

mysql> GRANT ALL PRIVILEGES ON userManagement.* to 'userManagement'@'localhost';
Query OK, 0 rows affected (0.00 sec)

#Install tomcat 
#tomcat7 requires java7
%sudo apt-get update
%sudo apt-get install tomcat7
#tomcat files in
 /var/lib/tomcat7 & /usr/share/tomcat7

#Configure tomcat
1. copy mysql-connector-java-5.1.34-bin.jar to /usr/share/tomcat/lib
2. configure tomcat basic auth to secure /user/admin/register web service:
add lines below to  /var/lib/tomcat7/conf/tomcat-users.xml
    <role rolename="admin"/>
    <user username="userManagementadmin" password="userManagementWSadmin" roles="admin"/>
3. inactive timeout in  /var/lib/tomcat7/conf/web.xml (default is 30mins)
    <session-config>
        <session-timeout>30</session-timeout>
    </session-config>
4. add the line below to tomcat /usr/share/tomcat7/bin/catalina.sh to increase PermSize:
JAVA_OPTS="$JAVA_OPTS -Xms256m -Xmx1024m -XX:PermSize=256m -XX:MaxPermSize=512m"
5. restart tomcat server
%sudo service tomcat7 restart

#Deploy userManagementWS to tomcat server
copy userManagementWS.war to /var/lib/tomcat7/webapps
the tomcat server will decompress the userManagementWS.war
#userManagementWS configurations
copy userManagementWS.xml to /var/lib/tomcat7/conf/Catalina/localhost/

#Deploy userManagement to tomcat server
copy userManagement.war to /var/lib/tomcat7/webapps
the tomcat server will decompress the userManagement.war

#have to change the url in /var/lib/tomcat7/webapps/userManagement/js/form.js  to        
var url = "http://edgeset.lanl.gov:8080/userManagementWS/user/uniqueEmail";

#change settings in /var/lib/tomcat7/webapps/userManagement/WEB-INF/classes/sys.properties

#setup admin user
#run script createAdminAccount.pl to add admin account with encrypted password to database
%perl createAdminAccount.pl -e admin@my.com -p admin

