<?php
require_once("PluploadHandler.php");

//read config
function read_config($configFile){
        $file_handle = fopen($configFile, "rb");
        while (!feof($file_handle) ) {
                $line_of_text = fgets($file_handle);
                $line_of_text=rtrim($line_of_text);
                if (preg_match("/=/",$line_of_text)){
                        $parts = explode('=', $line_of_text);
                        $array["$parts[0]"] = $parts[1];
                }
        }
        fclose($file_handle);
        return $array;
}

$edge_config=read_config(__DIR__."/../sys.properties");

if (!empty($edge_config["user_management"])){
        $session_return="false";
        $sid="";
        if (isset($_REQUEST["sid"])){
                $sid=$_REQUEST["sid"];
                //$session_return=shell_exec("ls /tmp/*$sid");
        }else{
                die('{"jsonrpc" : "2.0", "error" : {"code": 401, "message": "Invalid Session."}, "id" : "id"}');
                exit("Session Invalid");
        }
        //if (!empty($session_return) && ! preg_match("/$sid/i", $session_return)){
        //      die('{"jsonrpc" : "2.0", "error" : {"code": 102, "message": "Invalid Session."}, "id" : "id"}');
        //      exit("Session Invalid");
        //}
}

if (empty($_REQUEST["targetDir"]) || $_REQUEST["targetDir"] == "/" ) {
        exit("No target");
}

$domain=explode(".",$_SERVER['HTTP_HOST']);
$targetDir=$edge_config["edgeui_input"]."/$domain[0]";
if (file_exists("$targetDir")) {
        $targetDir = $edge_config["edgeui_input"]."/$domain[0]".$_REQUEST["targetDir"];
}else{
        $targetDir = $edge_config["edgeui_input"].$_REQUEST["targetDir"];
}

$maxDay = ($edge_config["edgeui_proj_store_days"]>0)? $edge_config["edgeui_proj_store_days"] : 1095;
$fileExt = ($edge_config["user_upload_fileext"])? $edge_config["user_upload_fileext"] : "fastq,fq,fa,fasta,fna,contigs,gbk,gbff,genbank,gb";
$maxFileAge = $maxDay * 24 * 60 * 60; // Temp file age in maxday days 

$ph = new PluploadHandler(array(
	'target_dir' => $targetDir,
	'allow_extensions' => $fileExt,
	'max_file_age' => $maxFileAge,
	'cb_check_file' => true
));
$ph->sendNoCacheHeaders();
$ph->sendCORSHeaders();
if ($result = $ph->handleUpload()) {
	die(json_encode(array(
		'OK' => 1,
		'info' => $result
	)));
} else {
	error_log($ph->getErrorMessage());
	die(json_encode(array(
		'OK' => 0,
		'error' => array(
			'code' => $ph->getErrorCode(),
			'message' => $ph->getErrorMessage()
		)
	)));
}

