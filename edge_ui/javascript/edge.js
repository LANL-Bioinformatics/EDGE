$( document ).ready(function()
{
	var debug=true;
	var pipeline="EDGE";
	var running_proj=0;
	var finished_proj=0;
	var failed_proj=0;
	var unstart_proj=0;
	var focusProjName; //id or name
	var focusProjRealName; // name
	var focusProjLogFile;
	var focusProjStatus;
	var focusProjType;
	var focusProjTime;
	var focusProjRunningStep;

	var interval = 1000*5; //check every 7 secs
	var updateProjInterval;
	var updateLogInterval;
	var inputFileID;
	var inputFileDir  = "/public/";
	//var username;
	//var password;
	var userType;
	var umSystemStatus = (localStorage.umStatus == 'true')? true : false;
	var umSystemURL = localStorage.umURL;

	var page = $( this );
	var allMainPage = $(".edge-main-page");
	
	//init page
	//$("#edge-content-upload-li").hide();
	$("#edge-phylo-ref-select-ref-div").hide();
	if (! umSystemStatus ){
		check_user_management();
	}else{
		$(".no-show-logout").hide();
		replaceUMurl();
		$('#edge-project-page-li').text('Public Project List');
		check_login();
	}
	
	$('#edge-user-btn').on("click",function(){
		$('#signInForm').popup('open');
	});

	edge_ui_init();
	sync_input(); //sync input with switch
	integrityCheck();
	updateProject(focusProjName);
	
	allMainPage.hide();
	$( "#edge-content-home" ).fadeIn();
	
	$( "a[href=#edge-content-home]" ).on( "click", function(){
		allMainPage.hide();
		$( "#edge-content-home" ).fadeIn();
		page.find( ".edge-navmenu-panel:not(.edge-panel-page-nav)" ).panel( "close" );
	});
	$( "a[href=#edge-content-pipeline]" ).on( "click", function(){
		allMainPage.hide();
		pipeline="EDGE";
		$("#edge-submit-info" ).children().remove();
		$("#edge-file-input-block").children().show();
		$(".btnAdd-edge-input").children().show();
		$(".edge-main-pipeline-input").show();
		$(".edge-qiime-pipeline-input").hide();
		$("#edge-content-pipeline" ).fadeIn("fast", function(){
			if (umSystemStatus && (localStorage.sid == "" || typeof localStorage.sid === "undefined") ){
				showWarning("Please login to run EDGE.");
			}
		});
		page.find( ".edge-navmenu-panel:not(.edge-panel-page-nav)" ).panel( "close" );
	});

	$( "a[href=#edge-qiime-pipeline]" ).on( "click", function(){
		allMainPage.hide();
		pipeline="qiime";
		$("#edge-submit-info" ).children().remove();
		$(".edge-main-pipeline-input").hide();
		$(".edge-qiime-pipeline-input").show();
		$("#edge-input-sequence").collapsible( "option", "collapsed", false );
		$('.edge-input-se-block').hide();
		$('.edge-input-pe-block').show();
		$('#edge-qiime-pipeline-dir-input').hide();
		$('#edge-qiime-rt-sw1').prop("checked",true).checkboxradio("refresh");
		$('#edge-qiime-rt-sw2').prop("checked",false).checkboxradio("refresh");
		$('#edge-qiime-rt-sw3').prop("checked",false).checkboxradio("refresh");
		integrityCheck();
		$("#edge-content-pipeline" ).fadeIn("fast", function(){
			if (umSystemStatus && (localStorage.sid == "" || typeof localStorage.sid === "undefined") ){
				showWarning("Please login to run EDGE.");
			}
		});
		page.find( ".edge-navmenu-panel:not(.edge-panel-page-nav)" ).panel( "close" );
	});

	$( "a[href=#edge-content-uploadfile]" ).on( "click", function(){
		allMainPage.hide();
		var maxFileSize = localStorage.maxFileSize || '10mb';
		$( "#edge-upload-maxFileSize" ).html(maxFileSize);
		$( "#edge-content-uploadfile" ).fadeIn("fast", function(){
			if (umSystemStatus && (localStorage.sid == "" || typeof localStorage.sid === "undefined") ){
				//destory file upload 
				var uploader = $("#uploader").pluploadQueue({
					buttons : {browse:false,start:false,stop:false}
				});
				$( "#edge-upload-maxFileSize").html("0 kb");
					showWarning("Please login to upload files.");
				}
        });
		page.find( ".edge-navmenu-panel:not(.edge-panel-page-nav)" ).panel( "close" );
    });

	
	$('#chck-rememberme').click(function() {
		if ($('#chck-rememberme').is(':checked')) {
             // save username and password	
            		localStorage.usrname = $('#signIn-email').val().replace(/ /g, '');
			localStorage.chkbx = $('#signIn-email').val();
		} else {
                	localStorage.usrname = '';
                	localStorage.chkbx = '';
        	}
        });
	 
	//if url
	function accessProjectFromUrl(qs) {
		qs = qs.split('+').join(' ');
		var params = {}, tokens,re = /[?&]?([^=]+)=([^&]*)/g;
		while (tokens = re.exec(qs)) {
			params[decodeURIComponent(tokens[1])] = decodeURIComponent(tokens[2]);
		}
		if( params.proj){
			updateProject(params.proj);
			updateReport(params.proj);
		}
	}
	accessProjectFromUrl(document.location.search);

	function check_user_management(){
    		$.ajax({
    			type: "POST",
			url: "./cgi-bin/edge_user_management.cgi", 
        		cache: false,
        		dataType: "json",
        	// send username and password as parameters to the Perl script
        		data: {"action": 'check','protocol': location.protocol, 'sid':localStorage.sid },
				error: function(XMLHttpRequest, textStatus, errorThrown) { 
					$('.no-show-logout').hide();
					$('#edge-projet-list-li').hide();
					$(".edge-user-btn").hide();
					$( "a[href=#edge-content-pipeline]" ).hide();
					$( "a[href=#edge-content-uploadfile]" ).hide();
					$( "a[href=#edge-qiime-pipeline]" ).hide();
					$('#edge-content-home').prepend("<h2 class='error'>Failed to check user management system. Please check server error log for detail or contact system administrator</h2>");
         			console.log("ERROR");
        		}, // error 
        		// script call was successful 
        		// data contains the JSON values returned by the cgi script 
        		success: function(data){
					localStorage.maxFileSize = data.maxFileSize;
					if (data.error) { // script returned error
						console.log(data.error);
						if (data.error.toLowerCase().indexOf("administrator") >= 0){
							$('.no-show-logout').hide();
							$('#edge-projet-list-li').hide();
							$(".edge-user-btn").hide();
							$( "a[href=#edge-content-pipeline]" ).hide();
							$( "a[href=#edge-content-uploadfile]" ).hide();
							$( "a[href=#edge-qiime-pipeline]" ).hide();
							$('#edge-content-home').prepend("<h2 class='error'>"+data.error+"</h2>")
						}else{
							// no configuration to use User management
							localStorage.umStatus = false;
							$('.no-show-logout').show();
							$(".edge-user-btn").hide();
							$("#action-unshare-btn").parent().hide();
							$("#action-share-btn").parent().hide();
							$("#action-publish-btn").parent().hide();
							uploadFiles(inputFileDir);
						}
          			} // if
          			else { // user_management is live
						umSystemStatus = true;
						localStorage.umStatus = umSystemStatus;
						localStorage.umURL = data.url;
						localStorage.sid = data.sid;
						umSystemURL = data.url;
						$(".edge-user-btn").show();;
						$(".no-show-logout").hide();
						$('.no-show-login').show();
						$('#edge-project-page-li').text('Public Project List');
						updateProject();
						replaceUMurl();
						if (data.socialLogin){
							localStorage.social = true;
							$('#edge-social-login').show();
							initSocial();
						}else{
							$('#edge-social-login').hide();
						}
						//	console.log(data.SUCCESS);
          				}
						check_login();
       				} // success
    		});
	}
	function replaceUMurl() {
		var umhref = (umSystemURL.indexOf(location.hostname)>0 || umSystemURL.indexOf("localhost")>0)? location.protocol  + "//" + location.host + "/userManagement": umSystemURL;
		$('#begin-password-reset').find('iframe').attr("src", umhref + "/resetPasswd.jsp");
		$('#signUpForm').find('iframe').attr("src", umhref + "/register.jsp");
		$('#UpdateProfileForm').find('iframe').attr("src", umhref + "/userUpdate.jsp");
	}
	function check_login(){
		if ( localStorage.sid && localStorage.sid !='' && localStorage.sid !='undefined' ){
    		$('.no-show-login').hide();
    		$(".no-show-logout").show();
			//username = localStorage.user;
			firstname = localStorage.fnname;
			//password = localStorage.pass;
			userType = localStorage.userType;
			FileTree("/" + localStorage.udir + "/");
			uploadFiles("/" + localStorage.udir + "/");

			$('#edge-user-btn').removeClass("ui-btn-icon-notext").addClass("ui-btn-icon-left edge-user-btn-login");
			$('#edge-user-btn').html(localStorage.fnname);
			$('#edge-project-page-li').text('My Project List');
			$('#edge-user-btn').unbind("click").on("click",function(){
				$('#popupUser').popup('open');
			});
		}else{
			FileTree(inputFileDir);
		}
	}

	//user menu
	$('#signInForm').popup(
		{afteropen: function(){
			if (localStorage.chkbx && localStorage.chkbx != '') {
				$('#chck-rememberme').prop('checked',true);
				$('#signIn-email').val(localStorage.usrname);
				$('#signIn-password').empty();
				$('#signIn-password').focus();
			} else {
				$('#chck-rememberme').prop('checked',false);
				$('#signIn-email').focus();
			}
		}},
		{positionTo:'window'},
		{transition:'fade'}
	);
	$('#popupUser').popup({positionTo:'#edge-user-btn'},{transition:'slidedown'});

	$('#signInForm').keypress(function (e) {
		var key = e.which;
		if(key == 13){  // the enter key code
    			$('#signIn-submit-btn').click();
    			return false;  
  		}
	});
	
	if ( localStorage.social ){
		$('#edge-social-login').show();
		initSocial();
	}else{
		$('#edge-social-login').hide();
	}
	function initSocial(){
		$.getScript( "javascript/social.js", function( data, textStatus, jqxhr ) {
			//console.log( textStatus ); // Success
		});
	}
	$('#edge-social-login').children('a').each(function( event ) {
		var network = $(this).html();
		$(this).on('click',function(){
			$('#signInForm').popup('close');
			var networkScope = "email";
			if (network.toLowerCase().indexOf("windows") >= 0){
				networkScope = "wl.emails";
			}
			if (network.toLowerCase().indexOf("linkedin") >= 0){
				networkScope = "r_emailaddress";
			}
			setTimeout(function(){
				$.mobile.loading( "show", {
					text: "Connect to " + network +" ...",
					textVisible: 1,
					html: ""
				});
			},10);
			setTimeout(function(){
				var HelloNetwork = hello(network);
				HelloNetwork.login({scope:networkScope}).then(function() {
				//	console.log('You are signed in to network');
					HelloNetwork.api('/me').then(function(r) {
				//		console.log(r);
						r.network = network;
						socialLogin(r);
					});
					$.mobile.loading( "hide" );
				}, function(e) {
					$.mobile.loading( "hide" );
				//	console.log('Signin error: ' + e.error.message);
					setTimeout(function(){$('#signInForm').popup('open');},100);
				});
			},100);
		});
	});
		
	function socialLogin(response){
		$.ajax({
			type: "POST",
			url: "./cgi-bin/edge_user_management.cgi",
			dataType: "json",
			cache: false,
			data: $.param(response) + '&' + $.param({"action": "sociallogin",'protocol': location.protocol, 'sid': localStorage.sid}),
			success: function(data){
				if (data.error){
					showWarning("Failed to use "+ response.network +"  acoount to login in. Please check server error log for detail.");
				}else{
					data.username = data.email;
					login(data);
				}
			},
			error: function(XMLHttpRequest, textStatus, errorThrown) {
				showWarning("Failed to use "+ response.network +" acoount to login in. Please check server error log for detail.");
			}
		});
	}

	$('#signIn-submit-btn').on('click', function(e){
		e.preventDefault();
		var data = {
			username : $('#signIn-email').val().replace(/ /g, ''), // get username
        		password : $('#signIn-password').val() // get password
		}
        	if (data.username.length == 0) {$('#signIn-email').addClass("highlight"); }  	
        	if (data.password.length == 0) {$('#signIn-password').addClass("highlight"); }  	
        	if (data.username && data.password) { // values are not empty
        		$('#signInForm').popup('close');
			login(data);
    		} // if
		//	setTimeout( function() { $('#dlg-invalid-credentials').popup('open').css('width','240px'); }, 300 );
	});
	function login (input){
		$.ajax({
			type: "POST",
			url: "./cgi-bin/edge_user_management.cgi", 
        		dataType: "json",
        		cache: false,
        		// send username and password as parameters to the Perl script
        		data: $.param(input) + '&' + $.param({"action": 'login','protocol': location.protocol}),
        		// script call was *not* successful
        		error: function(XMLHttpRequest, textStatus, errorThrown) {
					showWarning("Failed to login in. Please check server error log for detail.");
        		}, // error 
        		// script call was successful 
        		// data contains the JSON values returned by the cgi script 
        		success: function(data){
				localStorage.maxFileSize=data.maxFileSize;
				if (data.status) { // login success
					allMainPage.hide();
 				        $( "#edge-content-home" ).fadeIn();
 	        			$('.no-show-login').hide();
           				$(".no-show-logout").fadeIn("fast");
           				// add session storage
           				//localStorage.user = username;
					localStorage.fnname = data.firstname;
					localStorage.lnname = data.lastname;
           				//localStorage.pass = password;
					localStorage.userType = data.type;
					localStorage.sid = data.SESSION;
					localStorage.udir = data.UserDir;
					userType = data.type;
					FileTree("/" + data.UserDir + "/");
					uploadFiles("/" + data.UserDir + "/");

					$('#edge-user-btn').removeClass("ui-btn-icon-notext").addClass("ui-btn-icon-left edge-user-btn-login");
					$('#edge-user-btn').html(localStorage.fnname);
					$('#edge-project-page-li').text('My Project List');
					$('#edge-content-home').find(".error").remove();
					$('#popupUser').removeClass('highlight');
					$('#edge-user-btn').unbind("click").bind("click",function(){
						$('#popupUser').popup('open');
					});
					updateProject(focusProjName);
					var LoginMsg = true; 
					if (localStorage.LoginMsg){ LoginMsg =  ( localStorage.LoginMsg  === "true") ;}
					var cleanData = (data.CleanData > 0);
					if (  LoginMsg && cleanData ){
						showWarning("The intermediate bam/sam/fastq/gz files in the projects directory will be deleted if they are older than "+ data.CleanData  + " days. <br/><p><input type='checkbox' data-role='none' id='chk-close-warning'>Don't show this again.</p>");
						$('#chk-close-warning').click(function(){
							if ($('#chk-close-warning').is(':checked')) {
								localStorage.LoginMsg = false;
							}else{
								localStorage.LoginMsg = true;
							}
						});
					}

          			} // if
          			else { // login was failed
					password="";
          				$('#dlg-invalid-credentials').popup('open');
          				$('#dlg-invalid-credentials').on( "popupafterclose", function(){
							$('#signInForm').popup('open');
					});
          			} //else
       			} // success
      		}); // ajax
	}
	$('#resetPasswd-link').on('click', function(){
		$('#signInForm').popup('close');
		setTimeout( function() { $('#begin-password-reset').popup('open'); }, 300 );
	});
	$('#signUpBtn').on('click', function(){
		$('#signInForm').popup('close');
		setTimeout( function() { $('#signUpForm').popup('open'); }, 300 );
	});
	$('#signUp-submit-btn').on('click', function(){
		$('#signUpForm').popup('close');
		setTimeout( function() { $('#dlg-sign-up-sent').popup('open').css('width','480px'); }, 300 );
	});
	$('#UpdateProfileBtn').on('click', function(){
		$('#popupUser').popup('close');
		setTimeout( function() { 
			$('#UpdateProfileForm').popup('open').css('width','240px'); 
			$('#UpdateProfile-fn').val(localStorage.fnname);
			$('#UpdateProfile-ln').val(localStorage.lnname);
		}, 300 );
	});
	$('#updateProfile-cancel-btn').on('click',function(){
		$('#UpdateProfileForm').popup('close');
	});
	$('#updateProfile-submit-btn').on('click', function(){
		var newData = { 
			newPass: $('#UpdateProfile-password').val(),
			newFn:  $('#UpdateProfile-fn').val(),
			newLn: $('#UpdateProfile-ln').val(),
			action: 'update',
			protocol: location.protocol,
			sid: localStorage.sid
		};
		if($('#UpdateProfile-password').val() != $('#UpdateProfile-password-confirm').val()){
			$(".error").remove();
        		$('#UpdateProfile-password').addClass("highlight"); 
        		$('#UpdateProfile-password-confirm').addClass("highlight"); 
			$('#updateProfile-submit-btn').before("<p class='error'>Password not confirmed</p>")
		}else{
			$(".error").remove();
			updateProfile(newData);
		}
	});
	function updateProfile(input){
		$.ajax({
			type: "POST",
			url: "./cgi-bin/edge_user_management.cgi", 
        		dataType: "json",
        		cache: false,
        		data: input,
        		// script call was *not* successful
        		error: function(XMLHttpRequest, textStatus, errorThrown) {
				showWarning("Failed to udpate profile. Please check server error log for detail.");
        		}, // error 
        		success: function(data){
				if (data.error){
					showWarning("Failed to udpate profile. Please check server error log for detail.");
				}else{
					localStorage.fnname=input.newFn;
					localStorage.lnname=input.newLn;
					$('#UpdateProfileForm').popup('close');
					$('#edge-user-btn').html(input.newFn);
		    			$( "#edge_integrity_dialog_header" ).text("Message");
					$( "#edge_integrity_dialog_content" ).text("Your profile has been updated.");
					setTimeout( function() { $( "#edge_integrity_dialog" ).popup('open'); }, 300 );
				}
			} // success
      		}); // ajax
	};
	
	$('#signOutBtn').on('click', function(){
		logout("You have been logged out successfully.");
	});
	
	function logout(success_msg){
		$.ajax({
			type: "POST",
			url: "./cgi-bin/edge_user_management.cgi", 
        	dataType: "json",
        	cache: false,
        	data: {"action": 'logout','protocol': location.protocol, 'sid':localStorage.sid },
        	// script call was *not* successful
        	error: function(XMLHttpRequest, textStatus, errorThrown) {
				showWarning("Failed to login in. Please check server error log for detail.");
        	}, // error 
        	complete: function(data){
		   		$('#popupUser').popup('close');
				//localStorage.clear();
				localStorage.sid="";
				userType='';
				$('#signIn-password').empty();
				allMainPage.hide();
		 	    $( "#edge-content-home" ).fadeIn();
				// remove session, update button
				$('.no-show-login').fadeIn("fast");
		    	$(".no-show-logout").hide();
				$('#edge-user-btn').addClass("ui-btn-icon-notext").removeClass("edge-user-btn-login");
				$('#edge-user-btn').html("User");
				$('#edge-project-page-li').text('Public Project List');
				setTimeout( function() {updateProject(focusProjName);},5);
		    	$( "#edge_integrity_dialog_header" ).text("Message");
				$( "#edge_integrity_dialog_content" ).text(success_msg);
				setTimeout( function() { $( "#edge_integrity_dialog" ).popup('open'); }, 300 );
				$('#edge-user-btn').unbind("click").on("click",function(){
					$('#signInForm').popup('open');
				});
				FileTree(inputFileDir);
                                //destory file upload 
                                var uploader = $("#uploader").pluploadQueue({
                                    buttons : {browse:false,start:false,stop:false}
                                });
			} // success
      	}); // ajax
	};
	
	$(".edge-popup-close").on('click', function(e){
		e.preventDefault();
		$(this).parent().popup('close');
		location.href = location.href.split('#')[0];
		//$.mobile.navigate('/edge_ui/');
	});

	$('.edge-logo').find('h5').html('@'+location.hostname);
	$('.edge-logo').find('a').on('click',function(e){
        	e.preventDefault();
        	location.href = location.href.split('#')[0];
	});
	
	// initalize tooltipster
	$('.tooltip').tooltipster({
		theme:'tooltipster-light',
		maxWidth: '480',
		interactive: true,
	});
	// update qc tooltip content
	$('#qc-q-tooltip').tooltipster(
		'content', $('<span>Trim both end with the Phred quality. In brief, it will find the position in the read where trimming will end (argmax) based on the following equation: <img src="images/FaQCs.png"/> <br/> where l is the read length and Qu is the user-defined quality threshold, and trimming ends after the summation of Qu - Qi becomes negative.</span>')
	);
	$('.edge-sponsor-logo').tooltipster(
		'content', $('<h2>Mission</h2><p>EDGE addresses the critical rate limiting step in Next Generation Sequencing (NGS)  of genomic data analysis and enables OCONUS laboratories and resource restricted sites to utilize NGS technology for DoD&#39s SENSE and SHAPE missions.</p><p>EDGE is built to be a highly adaptable bioinformatics platform that is capable of rapidly analyzing and interpreting genomic sequence data into actionable results. This capability significantly reduces  the need for human resources and costs associated with genomic analysis for BWA detection and Biosurveillance efforts.</p>')
	);
	$('.edge-sponsor-logo').tooltipster('option','arrow',false);
	$("#edge-phylo-sra-acc-tooltip").tooltipster(
		'content', $('<span>(Internet requried) Input SRA accessions (comma separate for > 1 input) support studies (SRP*/ERP*/DRP*), experiments (SRX*/ERX*/DRX*), samples (SRS*/ERS*/DRS*), runs (SRR*/ERR*/DRR*), or submissions (SRA*/ERA*/DRA*).</span>')
	);
	$('.edge-aligner-options').tooltipster(
		'content', $('<span>Click <a href="http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#usage" target="_blank">Bowtie2 options</a> or <a href="http://bio-bwa.sourceforge.net/bwa.shtml#3" target="_blank">BWA mem options</a> for detail.</span>')
	);
	$('#edge-megahit-preset-l').tooltipster(
		'content', $('<span><table border="1"><tr><th>Presets</th><th>Targeting applications</th></tr><tr><td>meta</td><td>General metagenome assembly, such as guts</td></tr><tr><td>meta-sensitive</td><td>More sensitive metagenome assembly, but slower</td></tr><tr><td>meta-large</td><td>Large and complex metagenome assembly, such as soil</td></tr><tr><td>bulk</td><td>Experimental; assembly of standard bulk sequencing with sufficient depth</td></tr><tr><td>single-cell</td><td>Experimental; single-cell sequence assembly</td></tr></table></span>')
	);
	$("#edge-qiime-pro-tooltip").tooltipster(
		'content',$('<span><a href="images/pe_orientation.png" target="_blank"><img src="images/pe_orientation.png" width="460px"></img></a></span>')
	);
	
	$("#edge-qiime-mapping-file-tooltip").tooltipster(
		'content', $('<span>Metadata mapping files are used through-out QIIME, and provide per-sample metadata. The header for this mapping file starts with a pound (#) character, and generally requires a "SampleID", "BarcodeSequence", and a "Description", all tab separated. <a href="http://qiime.org/documentation/file_formats.html" target="_blank">Click here</a> for detail.</span>')
	);
	$("#edge-qiime-barcode-fq-file-tooltip").tooltipster(
		'content', $('<span>If the barcode has been removed from reads fastq file (demultiplexd.), please provide corresponding fastq file containing the barcode reads for each amplicon sequence. <a href="http://qiime.org/tutorials/processing_illumina_data.html" target="_blank">Click here</a> for detail.</span>')
	);
	$("#edge-qiime-UPARSE-tooltip").tooltipster(
		'content',$('<span>Use UPARSE pipeline clusters NGS amplicon reads into OTUs. <a href="http://drive5.com/uparse/" target="_blank">Edgar, R.C. (2013) Nature Methods</a></span>')
	);
	

	// filter serach callBack function
	function revealOrSearch ( index, searchValue ) {
		var ret = true;
		var idx;
		if (searchValue && searchValue.length > 2){  
			searchValue = searchValue.split(/[ ,]+/);
			var text = $(this).text().toLowerCase();
			var filttext = $(this).data("filtertext") || '';
			filttext = filttext.toLowerCase();
			for ( idx = 0 ; idx < searchValue.length ; idx++ ) {
				if (searchValue[idx] && searchValue[idx].length > 2){
					if( text.indexOf(searchValue[idx].toLowerCase(), 0) >= 0 || filttext.indexOf(searchValue[idx].toLowerCase(), 0) >= 0){
						ret = false; //not filter this one out
					}
	    	        	}
	    	    	}
	    	} 
	    	return ret;
	}
	function OrSearch( index, searchValue ) {
		var ret = false;
	    	var idx;
	    	var notMatch_cnt=0;
	    	var search_cnt=0;
		if (searchValue && searchValue.length > 0){  
			searchValue = searchValue.split(/[ ,]+/);
	    	    	var text = $(this).text().toLowerCase();
	    	    	var filttext = $(this).data("filtertext") || '';
	    	    	filttext = filttext.toLowerCase();
			for ( idx = 0 ; idx < searchValue.length ; idx++ ) {
				if (searchValue[idx] && searchValue[idx].length > 0){
					search_cnt++;
					if( text.indexOf(searchValue[idx].toLowerCase(), 0) < 0 && filttext.indexOf(searchValue[idx].toLowerCase(), 0) < 0){
						notMatch_cnt++;
					}
				}
	    	    	}
	    	    	if (search_cnt === notMatch_cnt) {
				ret = true;
	    	    	}
	    	} 
	    	return ret;
	}

	//options toggle
	$( ".edge-additional-options" ).hide();
	$( ".edge-additional-options-toggle" ).on( "click", function(){
		var id = $(this).attr('id');
		var optionsID = id.replace('toggle','options');
		$( "#" + optionsID ).toggle();
	});
	//ncbi sra input
	$( "#edge-sra-input-block" ).hide();
	$( ":radio[name='edge-sra-sw']" ).on("change",function(){
		if ( $(this).val() == 1 ){
			$( "#edge-sra-input-block" ).fadeIn('fast');
			$( "#edge-file-input-block" ).hide();
			$( ".btnAdd-edge-input" ).hide();
		}
		else{
			$( "#edge-sra-input-block" ).fadeOut('fast');
			$( "#edge-file-input-block" ).fadeIn('fast');
			$( ".btnAdd-edge-input" ).fadeIn('fast');
			$( "#edge-sra-acc" ).val('');
		}
	});
	//
	/*
    	// batch input 
    	$('#edge-batch-sample-input').click( function(e) {
    		e.preventDefault();
		var path = (umSystemStatus)? 'PublicData':'data';
		
    		var sampleInput = "#each unique project name in the bracket []\n" + "[Project1]\n" + "#q1=path/to/paired_end_file_1\n" + "q1=" + path + "/testData/Ecoli_10x.1.fastq\n" + 
    	                  "#q2=path/to/paired_end_file_2\n" + "q2=" + path + "/testData/Ecoli_10x.2.fastq\n" + "description=\"test batch input project 1\"\n";
    		sampleInput = sampleInput + "[Project2]\n" + "#s=path/to/single_end_file\n" + "s=" + path + "/testData/Ecoli_10x.1.fastq\n" + "description=\"test batch input project 2\"\n";
    		$('#edge-batch-text-input').val(sampleInput);
    		$('#edge-batch-text-input').textinput( "refresh" );
    		//$('#edge-batch-sample-input').hide();
    	});
    
	$('input[id=edge-batch-file-input]').on('change', readSingleFile);
 
	// Grab the text files and set them to our variable
	function readSingleFile(evt) {
    		//Retrieve the first (and only!) File from the FileList object
    		var f = evt.target.files[0]; 
		if (!f) {
        		alert("Failed to load file");
   		} else if (!f.type.match('text.*')) {
			showWarning(f.name + " is not a valid text file.");
		} else if (f.size > 5120) {
	         	// > 5Kb  It seems too big to process on the text 
			showWarning(f.name + " with size " + f.size + " bytes is too big (>5k)");
    	} else {
      		var r = new FileReader();
      		r.onload = function(e) { 
	    	var contents = e.target.result;
			$('#edge-batch-text-input').val(contents);
            	$('#edge-batch-text-input').textinput( "refresh" );
     		}
      		r.readAsText(f);
		} 
  	}
  	*/
  	
	//init reserch usage bar
	$( "input[id$='usage-bar']" ).each(function( event ) {
		$(this).parent().find("input").hide();
		$(this).parent().find(".ui-slider-handle").remove();
		$(this).parent().find(".ui-slider-track").css('margin','0 15px 0 15px').css('pointer-events','none');
	});

	//init file tree
	function FileTree(Dir){
		$( "#edge_file_tree" ).fileTree({
				root: Dir,
				script: './cgi-bin/jqueryFileTree.cgi',
			},function(file) {
				if ( ! /-dir-/.test(inputFileID)){
					var file_relative=file.replace(Dir, "");
					$( "#"+inputFileID ).val(file_relative);
					$( "#edge_file_dialog" ).popup('close');
				}
			},function(dire){
				//double click dir/files return dir
				if ( /-dir-/.test(inputFileID)){
					var file_relative=dire.replace(Dir, "");
					$( "#"+inputFileID ).val(file_relative);
					$( "#edge_file_dialog" ).popup('close');
				}
			}
		);
 
		//show file browser
		$( ".edge-file-selector" ).on( "click", function() {
			inputFileID = $(this).prevAll().children().prop("id");
		});
		$( "#edge_file_dialog" ).popup({
			beforeposition: function(){
				if ( /-dir-/.test(inputFileID)){
					$("#edge_file_dialog h4").text("Select a Directory");
				}else{
					$("#edge_file_dialog h4").text("Select a File");
				}
			}
		});
	}
	
	// bind event to content
	$( ".ui-btn.ui-input-clear" ).on( "click", function() {
		var dvalue = $(this).prev("input").prop( 'defaultValue' );
		$(this).prev("input").val(dvalue);
	});

	//button for adding input fileds
    	$('#btnAdd-edge-input-se, #btnAdd-edge-qiime-mapping-file').click( function(e) {
		e.preventDefault();
		var blockClass, blockIDprefix, label, inputID;
		var btnID = this.id;
		var delbtnClass = ".btnDel-edge-input";
		if ( btnID ==="btnAdd-edge-input-se") { 
			blockClass = ".edge-input-se-block";
			blockIDprefix = "edge-input-se-block";
			label = "Single-end FASTQ file";
			inputID = "edge-input-se";
		}
		else if ( btnID === "btnAdd-edge-qiime-mapping-file") { 
			blockClass = ".edge-qiime-pipeline-input-block";
			blockIDprefix = "edge-qiime-pipeline-input-block";
			label = "Mapping File";
			inputID = "edge-qiime-mapping-file-input";
		}
        	// how many "duplicatable" input fields we currently have
        	var num = $(blockClass).length;	
        	// the numeric ID of the new input field being added	
        	var newNum	= new Number(num + 1);		
        	var newElem = $('#' + blockIDprefix + num ).clone().attr('id', blockIDprefix + newNum);
        	newElem.find('label').attr( 'for', inputID + newNum ).text(label + '(' + newNum + ')');            
        	newElem.find('input').attr( 'id', inputID + newNum ).attr('name', inputID + '[]');
        	newElem.find(delbtnClass).css("visibility","visible");
        	// insert newElem
        	$('#' + blockIDprefix + num).after(newElem);
        	
        	// bind the selector 
        	newElem.find(".edge-file-selector").on( "click", function() {
			inputFileID = inputID + newNum;
		});
		    
		newElem.find(delbtnClass).on( "click", function() {
			$('#' + blockIDprefix + newNum ).remove();
			$('#' + btnID ).removeClass('ui-disabled');
		});
        	// business rule: limit the number of fields to 5
		if (newNum == 5) {
			$('#' + btnID ).addClass('ui-disabled');
			//alert('maximum fields reached')
        	}                        
    	});
    
   	$('#btnAdd-edge-input-pe').click( function(e) {
        	 
        	e.preventDefault();
        	// how many "duplicatable" input fields we currently have
        	var num = $('.edge-input-pe-block').length;	

        	// the numeric ID of the new input field being added	
        	var newNum	= new Number(num + 1);		
        	var newElem = $('#edge-input-pe-block' + num ).clone().attr('id', 'edge-input-pe-block' + newNum);
        	newElem.find('label:first').attr( 'for', 'edge-input-pe1-' + newNum ).text('Pair-1 FASTQ file (' + newNum + ')');                  
        	newElem.find('input:first').attr( 'id', 'edge-input-pe1-' + newNum ).attr('name', 'edge-input-pe1[]');
        	newElem.find('label:last').attr( 'for', 'edge-input-pe2-' + newNum ).text('Pair-2 FASTQ file (' + newNum + ')');    
        	newElem.find('input:last').attr( 'id', 'edge-input-pe2-' + newNum ).attr('name', 'edge-input-pe2[]');
        	newElem.find('.btnDel-edge-input-pe').css("visibility","visible");
        	// insert newElem
        	$('#edge-input-pe-block' + num).after(newElem);
        	
        	// bind the selector 
        	newElem.find(".edge-file-selector").on( "click", function() {
				inputFileID = $(this).prevAll().children().prop("id");
		    });
		    
		    newElem.find(".btnDel-edge-input-pe").on( "click", function() {
				$('#edge-input-pe-block' + newNum ).remove();
				$('#btnAdd-edge-input-pe' ).removeClass('ui-disabled');
		    });
		    
        	// business rule: limit the number of fields to 5
        	if (newNum == 5) {
        	     $('#btnAdd-edge-input-pe' ).addClass('ui-disabled');
        	     //alert('maximum fields reached')
        	}                        
    	});
  	$('.btnAdd-edge-phylo-ref-file, .btnAdd-edge-qiime-barcode-fq-file, .btnAdd-edge-ref-file').click( function(e) {
		e.preventDefault();
		var blockClass, blockIDprefix, inputID, label, selectClass;
		if ($(this).hasClass("btnAdd-edge-phylo-ref-file")){
			selectClass = ".btnAdd-edge-phylo-ref-file";
			blockClass = ".edge-phylo-ref-file-block";
			blockIDprefix = "edge-phylo-ref-file-block";
			inputID= "edge-phylo-ref-file";
			label = "Add Genome(s)";
		}else if ($(this).hasClass("btnAdd-edge-qiime-barcode-fq-file")){
			selectClass = ".btnAdd-edge-qiime-barcode-fq-file";
			blockClass = ".edge-qiime-barcode-fq-file-block";
			blockIDprefix = "edge-qiime-barcode-fq-file-block";
			inputID = "edge-qiime-barcode-fq-file-input";
			label = "Barcode Fastq File";
		}else if ($(this).hasClass("btnAdd-edge-ref-file")){
			selectClass = ".btnAdd-edge-ref-file";
			blockClass = ".edge-ref-file-block";
			blockIDprefix = "edge-ref-file-block";
			inputID = "edge-ref-file";
			label = "Reference Genome";
                }
		// how many "duplicatable" input fields we currently have
		var num = $(blockClass).length;
		console.log(num);	
		// the numeric ID of the new input field being added	
		var newNum	= new Number(num + 1);	
		var newElem = $('#' + blockIDprefix + num ).clone().attr('id', blockIDprefix + newNum);
		newElem.find('label:first').attr( 'for', inputID + newNum ).text(label + '(' + newNum + ')');
		newElem.find('input:first').attr( 'id', inputID + newNum ).attr('name', inputID);
		newElem.find(selectClass).remove();
		// insert newElem
		$('#' + blockIDprefix + num).after(newElem);
			// bind the selector 
			newElem.find(".edge-file-selector").on( "click", function() {
			inputFileID = $(this).prevAll().children().prop("id");
		});
		// business rule: limit the number of fields to 5
		if (newNum == 5) {
			$(selectClass).addClass('ui-disabled');
			//alert('maximum fields reached')
		}   
	});

	$( ".edge-collapsible-options > select" ).on( "change", function() {
		sync_input();
	});

	$( "#edge-all-on-btn" ).on( "click", function() {
		page.find(".edge-collapsible-options > select").val(1).slider("refresh");
		sync_input();
	});
	$( "#edge-all-exp-btn" ).on( "click", function() {
		page.find('div[data-role="collapsible"]').collapsible( "option", "collapsed", false );
	});
	$( "#edge-all-close-btn" ).on( "click", function() {
		page.find('div[data-role="collapsible"]').collapsible( "option", "collapsed", true );
	});
	
	$( "#edge-form-reset" ).on( "click", function() {
		page.find("form")[0].reset();
		 $('.ui-select select').val('').selectmenu('refresh');
		sync_input();
	});

	$( ".edge-collapsible-options div.ui-slider-switch" ).on( "mouseover", function() {
		$(this).parents('div[data-role="collapsible"]').collapsible("disable");
	});
	$( ".edge-collapsible-options div.ui-slider-switch" ).on( "mouseout", function() {
		$(this).parents('div[data-role="collapsible"]').collapsible("enable");
	});

	$( "#edge-view-log-btn" ).on( "mouseover", function() {
		getLog("."+focusProjLogFile);
	});

	$( "#edge_log_dialog" ).popup({
		beforeposition: function( event, ui ) {
			getLog("."+focusProjLogFile);
		},
		afteropen: function( event, ui ) {
			if( focusProjStatus != "finished" ){
				updateLogInterval = setInterval(
					function(){
						getLog("."+focusProjLogFile);
						if( focusProjStatus == "finished" ){
							clearInterval(updateLogInterval);
						}
					}, interval);
			}
		},
		afterclose: function( event, ui ) {
			clearInterval(updateLogInterval);
		}
	});

	//confirm dialog popup
	$( "a[id^='action']" ).on('click', function(){

		if (!focusProjName){
			showWarning("No focus project for the action");
			return false;
		}
		var action = $(this).attr("data");
		var actionContent = "Do you want to <span id='action_type'>"+action.toUpperCase()+"</span> project "+focusProjRealName+"?";

		//sample metadata
		if(action == "metadata-delete") {
			actionContent = "Do you want to <span id='action_type'>DELETE</span> project "+focusProjRealName+" sample metadata?";
		}
		if(action == "metadata-bsvedelete") {
			actionContent = "Do you want to <span id='action_type'>DELETE</span> project "+focusProjRealName+" sample metadata from BSVE?";
		}
		if(action == "metadata-bsveadd") {
			actionContent = "Do you want to <span id='action_type'>SHARE</span> project "+focusProjRealName+" sample metadata with BSVE?";
		}
		if(action == "metadata-bsveupdate") {
			actionContent = "Do you want to <span id='action_type'>UPDATE</span> project "+focusProjRealName+" sample metadata in BSVE?";
		}
		//END sample metadata


		if (action.indexOf("publish") < 0 && action.indexOf("bsve") < 0 ){
		//if (action.indexOf("publish") < 0){
			actionContent += "<p>This action can not be undone.</p>";
		}

		if (action == "share" || action == "unshare"){
			$("#edge_confirm_dialog_content").html( "<span id='action_type'>"+action.toUpperCase()+"</span> project " +focusProjRealName+ " to");
			setUserList(action,focusProjName);
		} 
		else if (action == "rename"){
				renameUserProject(action);
				//focusProjName is the projectID #
		}
		else{
			$("#edge_confirm_dialog_content").html(actionContent);
		}

		$("#edge_confirm_dialog a:contains('Confirm')").unbind('click').on("click",function(){
			actionConfirm(action,focusProjName);
		});
	});

	function renameUserProject(action){
		//send strings to #edge_cinfirm_diaglo_content to be printed in popup box
		var myVar = "Do you want to <span id='action_type'>"+action.toUpperCase()+"</span> "+focusProjRealName+"? <br> <br>";
		myVar +="Enter New Project Name <br>";
		myVar += "<input type= 'text' name='rename_project_Name' id='rename_project_Name_ID'> <br> <br>";
		myVar += "Enter New Project Description <br>";
		myVar += "<input type= 'text' name='rename_project_Desc' id='rename_project_Desc_ID'> <br> <br>";
		
		//sends to html page
		$("#edge_confirm_dialog_content").html(myVar);
	}

	function setUserList(action,pname){
		var actionContent= '<fieldset data-role="controlgroup" id="edge-userList" data-filter="true" data-filter-placeholder="Search users ...">';
		$.ajax({
			url:"./cgi-bin/edge_user_management.cgi",
			type: "POST",
			dataType: "json",
			cache: false,
			data: { "proj" : pname,  "action": action, 'protocol': location.protocol, 'sid':localStorage.sid },
			success: function(data){
				if (data.length>0){
					$("#edge_confirm_dialog a:contains('Confirm')").show();
					for (var i = 0 ; i < data.length ; i++){
						var FirstName = data[i].firstname;
						var LastName = data[i].lastname;
						var email = data[i].email;
						actionContent += "<input type='checkbox' data-mini='false' name='edge-userList' id='edge-userList_" + i  +"' value='" + email +  "'><label for='edge-userList_" + i  + "'>" + FirstName + " " + LastName + "</label>";
					}
					actionContent += '</fieldset>';
					$("#edge_confirm_dialog_content").append(actionContent);
					$('#edge-userList').filterable({
						children: ".ui-checkbox, label, input",
						filterCallback: revealOrSearch,
						create: function( event, ui ) {
   	         			 		$('#edge-userList').before("<div id='edge-userList-show-allnone-btn' style='margin-left:5%'> Show <a href='#' id='edge-userList-show-all'>All</a> | <a href='#' id='edge-userList-show-none'>None</a> </div>");
						},
						filter: function( event, ui ) {
							$('#edge-userList .ui-checkbox').children('label').each(function(){
								if($(this).hasClass('ui-checkbox-on')){
									$(this).removeClass('ui-screen-hidden');
									$(this).parent().removeClass('ui-screen-hidden');
									$('#edge-userList').controlgroup("refresh");
								}
							});
						}
					});
					$('#edge-userList-show-all').on('click',function(){
						$('#edge-userList .ui-checkbox').children().removeClass('ui-screen-hidden');
						$('#edge-userList .ui-checkbox').removeClass('ui-screen-hidden');
						$('#edge-userList').controlgroup("refresh");
					});
					$('#edge-userList-show-none').on('click',function(){
						$('#edge-userList .ui-checkbox').children().addClass('ui-screen-hidden');
						$('#edge-userList .ui-checkbox').addClass('ui-screen-hidden');
						$('#edge-userList').controlgroup("refresh");
					});
					
				}else{
					if (data.error){
						$("#edge_confirm_dialog_content").html(data.error);
					}else{
						$("#edge_confirm_dialog_content").html("All user have been "+ action + "d");
					}
				}
				$( "#edge_confirm_dialog" ).enhanceWithin().popup('open');
			},
			error: function(data){
				console.log(data);
				showWarning("ACTION FAILED: Please try again or contact your system administrator.");
			}
		});
	}

	$( "#edge_confirm_dialog" ).popup({
		afterclose: function( event, ui ) {
			updateProject(focusProjName);
		}
	});

	function actionConfirm (action,focusProjName,request){

		var rename_project = $("#rename_project_Name_ID").val(); //project name
		var project_description = $("#rename_project_Desc_ID").val(); //project desc


		var userChkArray=[];
		$('#edge-userList .ui-checkbox').children('label').each(function(){
			if($(this).hasClass('ui-checkbox-on')){
				userChkArray.push($(this).next().val());
			}
		});
		var shareEmail = userChkArray.join(',');
		//ajax sends the data to edge_action.cgi via Key, value
		var myAjaxRequest= $.ajax({
			url: "./cgi-bin/edge_action.cgi",
			type: "POST",
			dataType: "json",
			cache: false,
			data: { "proj" : focusProjName, "action": action, "userType":userType, "shareEmail" :shareEmail,'protocol': location.protocol, 'sid':localStorage.sid, "rename_project": rename_project, "project_description": project_description},
			beforeSend: function(){
				if (!request){
					$.mobile.loading( "show", {
						text: "Executing "+ action.toUpperCase() +" command...",
						textVisible: 1,
						html: ""
					});
				}
			},
			complete: function() {
				//$.mobile.loading( "hide" );
				//sample metadata
				if(action == 'metadata-delete') {
					updateReport(focusProjName);
				} 
				//NOTE: this change also affects other actions; need move it to if block?
				page.find( ".edge-action-panel" ).panel( "close" );
				//END sample metadata
			},
			success: function(data){
				$.mobile.loading( "hide" );
				if( data.STATUS == "SUCCESS" ){
					$( "#edge_integrity_dialog_header" ).text("Message");
					showWarning(data.INFO);
					if ( action == 'compare'){
						var loc = window.location.pathname;
						var edge_path = loc.substring(0,loc.lastIndexOf('/'));
						window.open(edge_path + data.PATH);
					}
				}
				else if( action == 'rename'){
					$( "#edge_integrity_dialog_header" ).text("Message");
					showWarning(data.INFO);
				} 
				else{
					showWarning(data.INFO);
				}
				if (action != 'compare'){
					updateProject(focusProjName);
				}
			},
			error: function(data){
				$.mobile.loading( "hide" );
				showWarning("ACTION FAILED: Please try again or contact your system administrator.");

			}
		});
		if (request){
			request.push(myAjaxRequest);
		}
	}

	//panel
	$( ".edge-navmenu-panel ul" ).listview();

	$( ".edge-navmenu-link" ).on( "click", function() {
		page.find( ".edge-navmenu-panel:not(.edge-panel-page-nav)" ).panel( "open" );
	});

	$( "div.edge-action-panel" ).panel({
 		beforeopen: function( event, ui ){
			updateProject(focusProjName);	
		}
	});

	$( ".edge-action-link" ).on( "click", function() {
		page.find( ".edge-action-panel" ).panel( "open" );
	});

	var w = $( "#edge-navmenu-panel" ).width();
	$( "#panel-folder" ).on( "click", function() {
		if( $( "#panel-folder-btn" ).hasClass("ui-icon-carat-l") ){
			$( "#panel-folder-btn" ).removeClass("ui-icon-carat-l").addClass("ui-icon-carat-r");
			$( "#edge-navmenu-panel" ).animate({ "left": "-="+w+"px" }, "fast" );
			$( "#edge-navmenu-panel" ).css("width","0px");
		}
		else{
			$( "#panel-folder-btn" ).removeClass("ui-icon-carat-r").addClass("ui-icon-carat-l");
			$( "#edge-navmenu-panel" ).animate({ "left": "+="+w+"px" }, "fast" );
			$( "#edge-navmenu-panel" ).css("width","25%");
		}
	});

	//update report
	var testKronaAnimation;
	function updateReport(pname) {
		$.ajax({
			url: "./cgi-bin/edge_report.cgi",
			type: "POST",
			dataType: "html",
			cache: false,
			data: { "proj" : pname, "sid":localStorage.sid , 'umSystem':umSystemStatus, 'protocol':location.protocol},
			beforeSend: function(){
				$.mobile.loading( "show", {
					text: "Load...",
					textVisible: 1,
					html: ""
				});
			},
			complete: function() {
				$.mobile.loading( "hide" );
			},
			success: function(data){
				//console.log("got response");
				allMainPage.hide();
				$( "#edge-content-report" ).html(data);
				if (  /The project does not exist/.test(data)){ 
					showWarning(data);
					return;
				}
				$( "#edge-content-report div[data-role='popup']" ).popup();
				$( "#edge-content-report > div[data-role='collapsible'] table " ).table();
				$( "#edge-content-report > div[data-role='collapsible']" ).collapsible();
				$( "#edge-content-report fieldset[data-role='controlgroup']" ).controlgroup();
				$( "#edge-content-report" ).show();
				$( "#edge-content-report" ).find("img").lazyLoadXT();
				$( "#edge-content-report" ).find("iframe").lazyLoadXT();
				$( "#edge-content-report" ).enhanceWithin();
				if (umSystemStatus && (localStorage.sid == "" || typeof localStorage.sid === "undefined")){
					$('#get_download_link').hide();
				}

				if ($( ".krona_plot" ).length){
                                        if(testKronaAnimation){clearInterval(testKronaAnimation);};
					testKronaAnimation = setInterval(function () {
						$( "#edge-content-report iframe" ).each(function(){
                                			var kf = this;
                                			var h = $("div[title='Help']", $(kf).contents());
                                			if ( $(h).size() ) {
                                        			$(h).parent().hide();
                                			}
                        			});
                			}, 100);
        			}else if(testKronaAnimation){
			
					clearInterval(testKronaAnimation);
				}
				
				$.getScript( "./javascript/edge-output.js" )
					.done(function( script, textStatus ) {
					//	console.log( "edge-output.js loaded: " + textStatus );
					})
					.fail(function( jqxhr, settings, exception ) {
						console.log( jqxhr, settings, exception );
					});
				var pstatus = $( "#edge-content-report" ).find("p:first").text().match( /Project Status: (.+)/m );

				if( pstatus[1] != "Complete" && pstatus[1] != "Archived"){
					$( "#edge-content-report div.ui-grid-a" ).hide();
					$( "#edge-content-report div.ui-grid-c" ).hide();
				}
				//console.log($( "#edge-get-contigs-by-taxa" ));
			},
			error: function(data){
				$.mobile.loading( "hide" );
				showWarning("Failed to retrieve the report. Please REFRESH the page and try again.");
			}
		});
	}

	//submit
	$( "#edge-form-submit" ).on( "click", function() {
		if (umSystemStatus && localStorage.sid == ""){
			showWarning("Please login to run EDGE.");
			return;
		}

		var projID=$("#edge-proj-name").val();
		$.ajax({
			url: "./cgi-bin/edge_submit.cgi",
			type: "POST",
			dataType: "json",
			cache: false,
			//data: $( "#edge-run-pipeline-form" ).serialize(),
			data: ( page.find("form").serialize() +'&'+ $.param({ 'pipeline': pipeline, 'protocol': location.protocol, 'sid':localStorage.sid })),
			beforeSend: function(){
				$(".list-info, .list-info-delete").fadeOut().remove(); //clear div
				page.find("input").removeClass("highlight");
				
				$.mobile.loading( "show", {
					text: "submitting...",
					textVisible: 1,
					html: ""
				});
			},
			complete: function(data) {
				$.mobile.loading( "hide" );
				$("#edge-submit-info").listview("refresh");
			},
			success: function(obj){
				// display general submission error
				$.each(obj, function(i,v){
					if( i!="PARAMS" && i!="SUBMISSION_STATUS" ){
						var dom;
						if( this.STATUS == "failure" ){
							dom = "<li data-icon='delete' data-theme='c' class='list-info-delete'><a href='#'>"+i+": "+this.NOTE+"</a></li>";
						}
						else{
							dom = "<li data-icon='info' class='list-info'><a href='#'>"+i+": "+this.NOTE+"</a></li>";
							if ( i == "PROJECT_NAME"){
								projID=this.NOTE.split(" ").pop();
							}
						}
						$( "#edge-submit-info" ).append(dom).fadeIn("fast");
					}
				});

				if( obj.SUBMISSION_STATUS == "success" ){
					$( ".highlight").removeClass("highlight");
					$( "#edge-submit-info" ).fadeIn("fast");
					var dom = "<li data-icon='info' class='list-info'><a href='#'>The job has been submitted successfully. Click to see open progress panel.</a></li>";
					
					$(dom).on( "click", function(){page.find( ".edge-action-panel" ).panel( "open" );})
					      .appendTo( "#edge-submit-info");
					updateProject(projID);
				}
				else{
					// display error information
					if(! $.isEmptyObject(obj.PARAMS)){
						$.each(obj.PARAMS, function(i,v){
							$("#"+i).addClass("highlight");
							$("#"+i).parents('div[data-role="collapsible"]').collapsible( "option", "collapsed", false );
							$( "#edge-submit-info" ).fadeIn("fast");
							var label = $("label[for='"+i+"']").text();
							var dom = "<li data-icon='delete' data-theme='c' class='list-info-delete'><a href='#'>"+label+": "+v+"</a></li>";
							
							$(dom).appendTo("#edge-submit-info").on("click", function(){
								$('html, body').animate({
									scrollTop: $("#"+i).offset().top-80
								}, 200);
							});
						});
					}
				}
				
				if( $(".list-info, .list-info-delete").size() ){
					var h = $( "#edge-submit-info" ).outerHeight();
					var h1 = parseInt(h)+100;
					$('html, body').animate({ scrollTop: "+="+h1+"px" }, 500);
				}

			},
			error: function(data){
				$( "#edge-submit-info" ).fadeIn("fast");
				var dom = "<li data-icon='delete' data-theme='c' class='list-info-delete'><a href='#'>FAILED to run submission CGI. Please check server error log for detail.</a></li>";
				$( "#edge-submit-info" ).append(dom).fadeIn("fast");
			}
		});
	});

	function integrityCheck() {
		// annotation and assembly
		$(":radio[name='edge-contig-taxa-sw']").on("change", function(){
			if ( $(this).val() == 1 ){
				if ( $( "#edge-assembly-sw" ).val() == 0 ){
					showWarning("Contig Classification function can only be performed when assembly function turns on.");
					$( "#edge-contig-taxa-sw2" ).prop('checked',true);
				}
			}
		});
		$(":radio[name='edge-primer-valid-sw']").on("change", function(){
			if ( $(this).val() == 1 ){
				if ( $( "#edge-assembly-sw" ).val() == 0 ){
					showWarning("Primer validation function can only be performed when assembly function turns on.");
					$( "#edge-primer-valid-sw2" ).prop('checked',true);
				}
			}
		});
		$(":radio[name='edge-primer-adj-sw']").on("change", function(){
			if ( $(this).val() == 1 ){
				if ( $( "#edge-assembly-sw" ).val() == 0 ){
					showWarning("Primer design function can only be performed when assembly function turns on.");
					$( "#edge-primer-adj-sw2" ).prop('checked',true);
				}
			}
		});
		$(":radio[name='edge-qc-sw']").on("change", function(){
			if ( $(this).val() == 0 ){
				if ( $("#edge-pp-sw").val() == 1 && $(":radio[name='edge-hostrm-sw']:checked").val()==0 ){
					showWarning("At least one function needs to be turned on!");
					$( "#edge-qc-sw1" ).click().checkboxradio("refresh");
				}
			}
		});
		$(":radio[name='edge-hostrm-sw']").on("change", function(){
			if ( $(this).val() == 0 ){
				if ( $("#edge-pp-sw").val() == 1 && $(":radio[name='edge-qc-sw']:checked").val()==0 ){
					showWarning("At least one function needs to be turned on!");
					$( "#edge-hostrm-sw1" ).click().checkboxradio("refresh");
				}
			}
		});

		$("#edge-assembly-sw").on("change", function(){
			if ( $(this).val() == 0 ){
				/*if ( $( "#edge-anno-sw1" ).is(':chekced') ){
					$( "#edge_integrity_dialog_content" ).text("Assembly function can not be turned off because annotation function is on.");
					$( "#edge_integrity_dialog" ).popup('open');
					$( "#edge-assembly-sw" ).val(1).slider("refresh");
				}
				elae*/ 
				if ( $( "#edge-primer-valid-sw1" ).is(':checked') ){
					$( "#edge-primer-valid-sw1" ).prop("checked",false).checkboxradio("refresh");
					$( "#edge-primer-valid-sw2" ).prop("checked",true).checkboxradio("refresh");
				}
				if ( $( "#edge-primer-adj-sw1" ).is(':checked') ){
					$("#edge-primer-adj-sw1").prop("checked",false).checkboxradio("refresh");
					$("#edge-primer-adj-sw2").prop("checked",true).checkboxradio("refresh");
				}
				if ( $( "#edge-contig-taxa-sw1" ).is(':checked') ){
					$("#edge-contig-taxa-sw1").prop("checked",false).checkboxradio("refresh");
					$("#edge-contig-taxa-sw2").prop("checked",true).checkboxradio("refresh");
				}
			}
		});
		$(":radio[name='edge-anno-sw']").on("change", function(){
			if ( $(this).val() == 0 ){
				if ( $( "#edge-orfs-sg-sw1" ).is(':checked') ){
					$( "#edge-orfs-sg-sw1" ).prop("checked",false).checkboxradio("refresh");
					$( "#edge-orfs-sg-sw2" ).prop("checked",true).checkboxradio("refresh");
				}
			}
		});
		$(":radio[name='edge-orfs-sg-sw']").on("change", function(){
			if ( $(this).val() == 1 ){
				if ( $( "#edge-anno-sw2" ).is(':checked') ){
					showWarning("ORFs Specialty Gene Profiling function can only be performed when annotation function turns on.");
					$( "#edge-orfs-sg-sw2" ).prop('checked',true);
				}
			}
		});
		$('#edge-anno-source').hide();
		
		$(":radio[name='edge-anno-tool']").on("change", function(){
			if($('#edge-anno-tool1').is(':checked')){
				$('#edge-anno-source').hide();
				$('#edge-anno-kingdom').show();
			}
			if($('#edge-anno-tool2').is(':checked')){
				$('#edge-anno-source').show();
				$('#edge-anno-kingdom').hide();
			}
		});

		$('#edge-spades-parameters').hide();
		$('#edge-megahit-parameters').hide();
		$(":radio[name='edge-assembler']").on("change", function(){
			if($('#edge-assembler1').is(':checked')){
				$('#edge-spades-parameters').hide();
				$('#edge-megahit-parameters').hide();
				$('#edge-idba-parameters').show();
			}
			if($('#edge-assembler2').is(':checked')){
				$('#edge-megahit-parameters').hide();
				$('#edge-idba-parameters').hide();
				$('#edge-spades-parameters').show();
			}
			if($('#edge-assembler3').is(':checked')){
				$('#edge-spades-parameters').hide();
				$('#edge-idba-parameters').hide();
				$('#edge-megahit-parameters').show();
			}
		});
		
		$('#edge-assembled-conti-file-div').hide();
		$(":radio[name='edge-assembled-contig-sw']").on("change", function(){
			if($('#edge-assembled-contig1').is(':checked')){
				$('#edge-assembler-sw').hide();
				$('#edge-idba-parameters').hide();
				$('#edge-spades-parameters').hide();
				$('#edge-assembled-conti-file-div').show();
			}
			if($('#edge-assembled-contig2').is(':checked')){
				$('#edge-assembler1').prop("checked",true).checkboxradio("refresh");
				$('#edge-assembler2').prop("checked",false).checkboxradio("refresh");
				$('#edge-assembler-sw').show();
				$('#edge-idba-parameters').show();
				$('#edge-assembled-conti-file-div').hide();
			}
		});
		
		$('.edge-input-se-block').hide();
		$('#edge-qiime-pipeline-dir-input').hide();
		$('#btnAdd-edge-input-se').hide();
		$(":radio[name='edge-qiime-rt-sw']").on("change", function(){
			if($('#edge-qiime-rt-sw1').is(':checked')){
				$('.edge-input-se-block').hide();
				$('.edge-input-pe-block').show();
				$('#edge-qiime-pro-sw').show();
				$('#btnAdd-edge-input-pe').show();
				$('#btnAdd-edge-input-se').hide();
				$('#edge-qiime-pipeline-dir-input').hide();
			}
			if($('#edge-qiime-rt-sw2').is(':checked')){
				$('.edge-input-pe-block').hide();
				$('.edge-input-se-block').show();
				$('#edge-qiime-pro-sw').hide();
				$('#btnAdd-edge-input-pe').hide();
				$('#btnAdd-edge-input-se').show();
				$('#edge-qiime-pipeline-dir-input').hide();
			}
			if($('#edge-qiime-rt-sw3').is(':checked')){
				$('.edge-input-se-block').hide();
				$('.edge-input-pe-block').hide();
				$('#edge-qiime-pipeline-dir-input').show();
				$('#btnAdd-edge-input-pe').hide();
				$('#btnAdd-edge-input-se').hide();
			}
		});
	}

	function sync_input() {	//disable inputs for default tools
		$( ".edge-collapsible-options > select" ).each( function() {
			var inputOpt = $(this).parents('div[data-role="collapsible"]');
			$(inputOpt).find("input").prop("disabled", false);
			$(inputOpt).find(".input-type-file select").prop("disabled", false);
			$(inputOpt).find(".input-type-file div").css("pointer-events", 'auto');
			
			if( $(this).val()==0 ){
				//$(inputOpt).find( "input:radio[value=1]" ).not("[name='edge-taxa-allreads']").prop("checked",false).checkboxradio( "refresh" );
				//$(inputOpt).find( "input:radio[value=0]" ).not("[name='edge-taxa-allreads']").prop("checked",true).checkboxradio( "refresh" );
				$(inputOpt).find( "input:radio").checkboxradio( "refresh" );
				$(inputOpt).find("input").prop("disabled", true);
				$(inputOpt).find(".input-type-file select").prop("disabled", true);
				$(inputOpt).find(".input-type-file div").css("pointer-events", 'none');
			}
		});
	};

	function getLog(texturl) {
		$.ajax({
			url: texturl,
			dataType: 'text',
			cache: false,
			success: function(text) {
				$("#edge-log-view").html(text);
				
				//auto scroll down if log-auto-scroll button is on
				if( $("#log-auto-scroll").is(':checked') ){
					var height = $("#edge_log_dialog pre").prop('scrollHeight');
					$("#edge_log_dialog pre").scrollTop(height);
				}
			},
			error: function(text) {
				$("#edge-log-view").text("No log retrieved from "+texturl);
			}
		})
	};

	function showWarning( dialog_content ) {
		$( "#edge_integrity_dialog_content" ).html(dialog_content);
		$( "#edge_integrity_dialog" ).popup('open');
	}

	$("#edge-project-page-li").on("click",function(){
		updateProjectsPage('user');
	});
	function updateProjectsPage(view) {
		$.ajax({
			url: "./cgi-bin/edge_projectspage.cgi",
			type: "POST",
			dataType: "html",
			cache: false,
			data: {'umSystem':umSystemStatus,'userType':userType,'view':view,'protocol': location.protocol, 'sid':localStorage.sid},
			beforeSend: function(){
				$.mobile.loading( "show", {
					text: "Load...",
					textVisible: 1,
					html: ""
				});
			},
			complete: function() {
				$.mobile.loading( "hide" );
			},
			success: function(data){
				updateProject();
				allMainPage.hide();
				$( "#edge-project-page" ).html(data);
				$( "#edge-project-page" ).show();
				$( ".edge-project-page-link").unbind('click').on('click', function(e){
					e.preventDefault();
					var pname = $(this).attr("data-pid");
					updateReport(pname);
					updateProject(pname);
				});
				$(".tooltip").tooltipster({
					theme:'tooltipster-light',
					maxWidth: '480',
					interactive: true,
					multiple:true
				});
				$("#edge-projpage-action").children('a').on("click", function(){
					var action = $(this).val();
					var actionContent = "Do you want to <span id='action_type'>"+action.toUpperCase()+"</span> projects " ;
					//sample metadata
					if(action === "metadata-bsveadd") {
						actionContent = "Do you want to <span id='action_type'>SHARE</span> projects' metadata with BSVE?" ;
					}
					//END sample metadata
					if ( action === "0" ){
						return;
					}
					if ( action  === "show-all" ){
						updateProjectsPage('admin');
						return;
					}else{
						if ( $('[name="edge-projpage-ckb"]:checked').length === 0 ){
							showWarning("There are no projects selected.");
							return;
						}
						if ( $('[name="edge-projpage-ckb"]:checked').length < 2 && action === "compare"){
							showWarning("There is only one project selected for comparison.");
							return;
						}
						var projnames=[];
						var projids=[];
						$('[name="edge-projpage-ckb"]:checked').each(function( event ) {
							projnames.push("<li>"+$(this).closest('td').next('td').find('.edge-project-page-link').text()+"</li>");
							if ( action === 'compare'){
								projids.push($(this).val());
							}else{
								projids.push($(this).closest('td').next('td').find('.edge-project-page-link').attr('data-pid'));
							}
						});
						actionContent += '<ul>' + projnames.join('\n') + '</ul>';
						var focusProjCodes = projids.join();
						$("#edge_confirm_dialog_content").html(actionContent);
						$( "#edge_confirm_dialog" ).enhanceWithin().popup('open');
						if ( action === "share" ){
							setUserList(action,focusProjCodes);
						}
						$("#edge_confirm_dialog a:contains('Confirm')").unbind('click').on("click",function(){
							if ( action === "compare" ){
								actionConfirm(action,focusProjCodes);
								//updateProjectsPage('user');
							}else{
								var actionRequest=[];
								for (var i = 0; i < projids.length; i++) { 
									$.mobile.loading( "show", {
										text: "Executing "+ action.toUpperCase()  + " command...",
										textVisible: 1,
										html: ""});
									actionConfirm(action,projids[i],actionRequest);
								}
								// wait for all ajax request done
								$.when.apply(null,actionRequest).done(function(){
									$( "#edge_integrity_dialog" ).popup('close');
									$( "#edge_integrity_dialog_header" ).text("Message");
									$( "#edge_integrity_dialog" ).popup("reposition",{positionTo: 'window'});
									//updateProjectsPage('user');
									//sample metadata
									if(action === "metadata-bsveadd") {
										showWarning('Metadata of projects' + '<ul>' + projnames.join('\n') + '</ul>'+ 'have been shared with BSVE.');
									} else {
									//END sample metadata
										showWarning('Projects' + '<ul>' + projnames.join('\n') + '</ul>'+ 'have been ' + action +'d');
									}
									//showWarning('Projects' + '<ul>' + projnames.join('\n') + '</ul>'+ 'have been ' + action +'d');
								});
							}
						});
					}
					
				});
				$("#edge-project-list-filter").filterable({
					children: $( "#edge-project-list-filter" ).find("tbody").children(),
					create: function (event,ui){
						 $('#edge-project-page').find('.ui-filterable').addClass('ui-grid-a').css('margin','0px');
						 $('#edge-project-page').find('.ui-filterable').find('div').addClass('ui-block-b');
						 $('#edge-projpage-action').addClass('ui-block-a');
						 $('#edge-project-page').find('.ui-block-b').before($('#edge-projpage-action'));
					},
					filterCallback: OrSearch,
					filter:function( event, ui ) { 
						$("#edge-project-list-filter .ui-collapsible").each(function(){
							if($(this).find('tbody').children().not(".ui-screen-hidden").size() >0){
								$(this).collapsible("option", "collapsed", false);		
							}else{
								$(this).collapsible("option", "collapsed", true);
							}
						});
					}
				});
				$( "#edge-project-page" ).enhanceWithin();

				$( "#edge_confirm_dialog" ).popup({
					afterclose: function( event, ui ) {
						updateProject(focusProjName);
					}
				});
			/*	//use tablesortr library 
				$("#edge-project-page-table").tablesorter(
					{
						widthFixed: true, 
						widgets: ['zebra','filter'],
						widgetOptions : {
						 	reflow_className    : 'ui-table-reflow',
						 	reflow_headerAttrib : 'data-name',
							reflow_dataAttrib   : 'data-title',
							filter_reset: '.reset'
						}
					}
				);*/
			},
			error: function(data){
				$.mobile.loading( "hide" );
				$( "#edge_integrity_dialog_content" ).text("Failed to retrieve the report. Please REFRESH the page and try again.");
				$( "#edge_integrity_dialog" ).popup('open');
			}
		});
	};


	function edge_ui_init () {
		$.ajax({
			url: './cgi-bin/edge_info.cgi',
			type: "POST",
			dataType: "json",
			cache: false,
			data: { "init" : '1', 'umSystem':umSystemStatus, 'protocol':location.protocol, 'sid':localStorage.sid },
			complete: function(data){
			},
			success: function(obj) {
				if( obj.INFO.UPLOAD == "false"){
					var uploader = $("#uploader").pluploadQueue({
					buttons : {browse:false,start:false,stop:false}
					});
					$("#edge-content-upload-li").hide();
				}else{
					$("#edge-content-upload-li").show();
				}
				( obj.INFO.ARCHIVE == "true" )?	$("#action-archive-btn").show():$("#action-archive-btn").hide();
				
				if ( obj.INFO.MQC == "true" ){
					$("#edge-pp-sw").closest('div[data-role="collapsible"]').show();
				}else{
					$("#edge-pp-sw").val(0).slider("refresh");
					$("#edge-pp-sw").closest('div[data-role="collapsible"]').hide();
				}
				if( obj.INFO.MAA == "true" ){
					$("#edge-assembly-sw").closest('div[data-role="collapsible"]').show();
				}else{
					$("#edge-assembly-sw").val(0).slider("refresh");
					$("#edge-assembly-sw").closest('div[data-role="collapsible"]').hide();
				}
				if( obj.INFO.MRBA == "true"){
					$("#edge-ref-sw").closest('div[data-role="collapsible"]').show();
				}else{
					$("#edge-ref-sw").val(0).slider("refresh");
					$("#edge-ref-sw").closest('div[data-role="collapsible"]').hide();
				}
				if( obj.INFO.MTC == "true" ){
					$("#edge-taxa-sw").closest('div[data-role="collapsible"]').show();
				}else{
					$("#edge-taxa-sw").val(0).slider("refresh");
					$("#edge-taxa-sw").closest('div[data-role="collapsible"]').hide();
				}
				if( obj.INFO.MPA == "true" ){
					$("#edge-phylo-sw").closest('div[data-role="collapsible"]').show();
				}else{
					$("#edge-phylo-sw").val(0).slider("refresh");
					$("#edge-phylo-sw").closest('div[data-role="collapsible"]').hide();
				}
				if( obj.INFO.MSGP == "true" ){
					$("#edge-sg-sw").closest('div[data-role="collapsible"]').show();
				}else{
					$("#edge-sg-sw").val(0).slider("refresh");
					$("#edge-sg-sw").closest('div[data-role="collapsible"]').hide();
				}
				if( obj.INFO.MPPA == "true"){
					$("#edge-primer-sw").closest('div[data-role="collapsible"]').show();
				}else{
					$("#edge-primer-sw").val(0).slider("refresh");
					$("#edge-primer-sw").closest('div[data-role="collapsible"]').hide();
				}
				(  obj.INFO.MQIIME == "true")?$( "a[href=#edge-qiime-pipeline]" ).show():$( "a[href=#edge-qiime-pipeline]" ).hide();
				
				if( String(obj.INFO.UMSYSTEM) != String(localStorage.umStatus) ){
					check_user_management();
				}
				$('#edge-proj-cpu').val(obj.INFO.RUNCPU);
			},
			error: function(data){
				showWarning("Failed to initialized EDGE GUI interface. Please check server error log for detail.");
			}
		});
	};

	function updateProject(pname) {
		focusProjName = pname;
		$.ajax({
			url: './cgi-bin/edge_info.cgi',
			type: "POST",
			dataType: "json",
			cache: false,
			data: { "proj" : pname, 'umSystem':umSystemStatus, 'protocol':location.protocol, 'sid':localStorage.sid },
			complete: function(data){
				/*
				console.log("finished_proj="+finished_proj);
				console.log("running_proj="+running_proj);
				console.log("failed_proj="+failed_proj);
				console.log("focusProjName="+focusProjName);
				console.log("focusProjStatus="+focusProjStatus);
				*/
				$( "#edge-project-list-ul > li" ).off('click').on("click", function(e){
					e.preventDefault();
					var pname = $(this).children("a").attr("data-pid");
					updateReport(pname);
					updateProject(pname);
				});
				clearInterval(updateProjInterval);
				updateProjInterval = setInterval(function(){ if (running_proj>0 || unstart_proj>0){updateProject(pname)}}, interval);
			},
			success: function(obj) {
				// update focused proj
				focusProjName   = obj.INFO.NAME;
				focusProjRealName   = obj.INFO.PROJNAME || obj.INFO.NAME;
				focusProjLogFile = obj.INFO.PROJLOG;
				focusProjStatus = obj.INFO.STATUS;
				focusProjTime   = obj.INFO.TIME;
				focusProjType	= obj.INFO.PROJTYPE;
				projListNumShow	= obj.INFO.PROJLISTNUM;
				finished_proj   = 0;
				running_proj    = 0;
				failed_proj     = 0;
				unstart_proj	= 0;
				//force logout
				if( obj.INFO.SESSION_STATUS === "invalid" && localStorage.sid ){
					logout("Session expired. Please login again.");
					return;
				}

				//resource usage bar
				$("#cpu-usage-bar").val(obj.INFO.CPUU).slider("refresh");
				$("#mem-usage-bar").val(obj.INFO.MEMU).slider("refresh");
				$("#disk-usage-bar").val(obj.INFO.DISKU).slider("refresh");
				$("#cpu-usage-val").html(obj.INFO.CPUU+" %");
				$("#mem-usage-val").html(obj.INFO.MEMU+" %");
				$("#disk-usage-val").html(obj.INFO.DISKU+" %");
				//console.log(focusProjType);
				if (focusProjType){
					(focusProjType.toLowerCase().indexOf("shared") >= 0)?
							$("#action-unshare-btn").parent().show():$("#action-unshare-btn").parent().hide();
					(focusProjType.toLowerCase().indexOf("publish") >= 0)?
							$("#action-publish-btn").attr('data','unpublish').text("Make project private"):
							$("#action-publish-btn").attr('data','publish').text("Make project public");
					(focusProjType.toLowerCase().indexOf("guest") >= 0 && !debug)?
							$("#action-share-btn").parent().parent().hide():
							$("#action-share-btn").parent().parent().show();
				}

				//sample metadata
				focusProjShowMeta = obj.INFO.SHOWMETA;
				focusProjIsOwner = obj.INFO.ISOWNER;
				focusProjHasMeta = obj.INFO.HASMETA;
				focusProjMetaBSVE = obj.INFO.METABSVE;

				if(focusProjShowMeta && (focusProjIsOwner || userType == 'admin')) {
					$("#action-metadata-bsve-submit-btn").parent().parent().show();
					if(focusProjHasMeta) {
							$("#action-metadata-delete-btn").show();
							$("#edge-pg-metadata-edit").attr('data','edit').text("Edit");
							$("#action-metadata-bsve-submit-btn").show();
					} else {
							$("#edge-pg-metadata-edit").attr('data','edit').text("Add");
							$("#action-metadata-delete-btn").hide();
							$("#action-metadata-bsve-delete-btn").hide();
							$("#action-metadata-bsve-submit-btn").hide();
					}
					(focusProjMetaBSVE)?
							$("#action-metadata-bsve-delete-btn").show():
							$("#action-metadata-bsve-delete-btn").hide();
					(focusProjMetaBSVE)?
							$("#action-metadata-bsve-submit-btn").attr('data','metadata-bsveupdate').text("Update in BSVE"):
							$("#action-metadata-bsve-submit-btn").attr('data','metadata-bsveadd').text("Share with BSVE");
				} else {
					$("#action-metadata-bsve-submit-btn").parent().parent().hide();
				}
				//END sample metadata

				// project list
				if(! $.isEmptyObject(obj.LIST)){
					$( "#edge-project-list-ul .edge-proj-list-li" ).remove();
		
					var listIdOrder = Object.keys(obj.LIST);
					listIdOrder.sort(function(a,b){ return obj.LIST[a].TIME < obj.LIST[b].TIME ? -1 : obj.LIST[a].TIME > obj.LIST[b].TIME; }).reverse();
					if (projListNumShow == 0){ projListNumShow = 9999999;}
					var proj_count=0;
					$.each(listIdOrder, function(i,v){
						var proj_list_obj = obj.LIST[v];
						if( proj_list_obj.NAME ){
							proj_count++;
							// with user management, NAME becomes unique project id
							var projname = proj_list_obj.PROJNAME;
							var name   = proj_list_obj.NAME;
							var time   =  proj_list_obj.TIME;
							var pstatus = proj_list_obj.STATUS;
							if (!projname) { projname = name;}
							var bgColor = "#4d4d4d";
							var displayList = "hiddenProjList"; 
							var desc = proj_list_obj.DESC || "No description";
							desc = desc + " ("+pstatus+")";
									
							switch ( pstatus ) {
								case "finished":
									projClass = "edge-time-bg-green";
									projIcon  = "ui-icon-check";
									finished_proj++;
									break;
								case "running":
									projClass = "edge-time-bg-orange";
									projIcon  = "ui-icon-load";
									running_proj++;
									break;
								case "archived":
									projClass = "edge-time-bg-black";
									projIcon = "ui-icon-arrow-u-r";
									finished_proj++;
									break;
  								case "failed":
									projClass = "edge-time-bg-red";
									projIcon  = "ui-icon-delete";
									failed_proj++;
									break;
  								default:
  									projClass = "edge-time-bg-grey";
  									projIcon  = "ui-icon-refresh";
									unstart_proj++;
							}
							if (focusProjName == name){
								bgColor = "#041B8E";
							}
							if (proj_count < projListNumShow){
								displayList = "";
							}
							var dom = "<li class='edge-proj-list-li "  + displayList + "'><div class='edge-project-time "+projClass+"'>"+time+"</div><a href='' style='background:"+bgColor+" !important' class='edge-project-list ui-btn ui-btn-icon-right "+projIcon+"' title='"+desc+"' data-pid='"+name+"'>"+projname+"</a></li>";
							$(dom).appendTo( "#edge-project-list-ul" );
						}
					});
				}else{
					$( "#edge-project-list-ul .edge-proj-list-li" ).remove();
				}
				if( $( ".edge-proj-list-li" ).size() == 0 ){
					var dom = "<li class='edge-proj-list-li ui-disabled'><a href='#' class='edge-project-list ui-btn ui-btn-icon-right ui-icon-check'>No project found</a></li>";
					$( "#edge-project-list-ul" ).append(dom);
				}

				$('.hiddenProjList').hide();
				$('#edge-project-list-ul').filterable({
					children:".edge-proj-list-li li",
					filterCallback:OrSearch,
					filter:function(event,ui){
						if ($('#filterProjectList-input').val()){
							$('.hiddenProjList').not('ui-screen-hidden').show();
						}
					}
				});

				// progress info
				if(! $.isEmptyObject(obj.PROG))
				{
					$( "#edge-progress-ul > li" ).fadeOut().remove();
					var projname = obj.INFO.PROJNAME || obj.INFO.NAME;
					var dom = "<li data-role='list-divider'>"+projname+"</li>";
					$( "#edge-progress-ul" ).append(dom);
		
					var listIdOrder = Object.keys(obj.PROG);
					listIdOrder.sort(function(a,b){ return parseInt(a) - parseInt(b);});

					focusProjRunningStep=100;

					$.each(listIdOrder, function(i,v){
						var prog_list_obj = obj.PROG[v];
						var name   = prog_list_obj.NAME;
						var pstatus = prog_list_obj.STATUS;

						switch ( pstatus ) {
							case "skip":
								projIcon = (parseInt(v) > focusProjRunningStep) ? "ui-icon-forbidden" : "ui-icon-forbidden edge-icon-bg-green";
								projText = (parseInt(v) > focusProjRunningStep) ? "Off" : "Skipped";
								break;
							case "finished":
								projIcon = "ui-icon-forward edge-icon-bg-green";
								projText = "Result exists. Skipped this step.";
								break;
							case "running":
								focusProjRunningStep = parseInt(v);
								projIcon = "ui-icon-recycle edge-icon-bg-orange";
								projText = "Running";
								break;
							case "done":
								projIcon = "ui-icon-check edge-icon-bg-green";
								projText = "Complete";
								break;
							case "failed":
								focusProjRunningStep = parseInt(v);
								projIcon = "ui-icon-delete edge-icon-bg-red";
								projText = "Failed";
								break;
							default:
								projIcon = "ui-icon-clock";
								projText = "Incomplete";
						}

						dom = "<li><a href='#' class='ui-btn ui-btn-icon-right "+projIcon+"' title='"+projText+"'>"+name+"</a></li>";
						if ( pstatus != "skip"){
							$( "#edge-progress-ul" ).append(dom);
						}
					});

					dom = "<li data-role='list-divider' class='edge-proj-last-check'>Last checked: "+obj.INFO.TIME+"</li>";

					//all incomplete
					if( focusProjRunningStep == 100 && focusProjStatus != "finished"){
						 $( "#edge-progress-ul" ).find(".ui-icon-forbidden").removeClass("edge-icon-bg-green");
					}

					$( "#edge-progress-ul" ).append(dom);
					$( "#edge-submit-info" ).listview("refresh");
					$( "#edge-progress-ul" ).listview("refresh");
					$( "#edge-project-list-ul" ).listview("refresh");
				}else{
					$( "#edge-progress-ul > li" ).remove();
				}
			},
			error: function(data){
				$( "#edge-submit-info" ).fadeIn("fast");
				var dom = "<li data-icon='delete' data-theme='c' class='list-info-delete'><a href='#'>FAILED to retrieve project info. Please check server error log for detail.</a></li>";
				$( "#edge-submit-info" ).append(dom);
			}
		});
	};





        //upload files  
        function uploadFiles(userDir){
		var target = (umSystemStatus)? "MyUploads/":"Uploads/";
		var maxFileSize = localStorage.maxFileSize || '10mb';
                var uploader = $("#uploader").pluploadQueue({
                    // General settings
                    runtimes : 'html5,flash,silverlight,html4',
                    url : './cgi-bin/upload.php?targetDir='+userDir+target,

                    // User can upload no more then 20 files in one go (sets multiple_queues to false)
                    max_file_count: 20,

                    chunk_size: '1mb',

                    // Resize images on clientside if we can
                    resize : {
                            width : 200,
                            height : 200,
                            quality : 90,
                            crop: true // crop to exact dimensions
                    },

                    filters : {
                            // Maximum file size
                            max_file_size :  maxFileSize,
                            // Specify what files to browse for
                            mime_types: [
                                    {title : "text/plain", extensions : "fastq,fq,fa,fasta,fna,contigs,gbk,genbank,gb,txt,text,config,ini,xls,xlsx"},
                                    {title : "application/x-gzip", extensions : "gz"},
                            ]
                    },

                    // Rename files by clicking on their titles
                    rename: true,

                    // Sort files
                    sortable: true,

                    // Enable ability to drag'n'drop files onto the widget (currently only HTML5 supports that)
                    dragdrop: true,

                    //prevent_duplicates: true,

                    // PreInit events, bound before any internal events
                    preinit: {
                       Init: function(up, file) {
                          if (umSystemStatus && (localStorage.sid == "" || typeof localStorage.sid === "undefined") ){
                             //destory file upload 
                             var uploader = $("#uploader").pluploadQueue({
                                 buttons : {browse:false,start:false,stop:false}
                             });
                           }
                       }
                    }, 

                    // Post init events, bound after the internal events
                    init : {
                       BeforeUpload: function (up, file) {
			  if (! umSystemStatus){
				return true;
			  }
                          var url = "./cgi-bin/plupload_session.cgi?sid="+localStorage.sid;
                          $.ajax({
                             type: "GET",
                             url: url,
                             dataType: "text",
                             contentType: "application/text; charset=utf-8",
                             cache: false
                          }).done(function (msg) {
                             msg = msg.replace(/\s+/g,'');
                             if (msg == 'true') {
                                file.status = plupload.UPLOADING;
                                up.trigger("UploadFile", file);
                             } else {
                                logout("Session expired. Please login again.");
                                return false;
                             }
                          });
                          return false;
                       },

                       UploadComplete: function(up, files) {
                          // Called when all files are either uploaded or failed
			  (umSystemStatus)?FileTree(userDir):FileTree(inputFileDir);
                       }
                    },

            });
    };

        // Handle the case when form was submitted before uploading has finished
        $('#edge-uploadfile-form').submit(function(e) {

                if (umSystemStatus && localStorage.sid == ""){
                        showWarning("Please login to upload files.");
                        return;
                }

                // Files in queue upload them first
                if ($('#uploader').plupload('getFiles').length > 0) {

                        // When all files are uploaded submit form
                        $('#uploader').on('complete', function() {
                                $('#edge-uploadfile-form')[0].submit();
                        });

                        $('#uploader').plupload('start');
                } else {
                        alert("You must have at least one file in the queue.");
                }
                return false; // Keep the form from submitting
        });

	// scroll top
	$(window).scroll(function() {
		if ($(this).scrollTop()>600) {
			$('#scroll-up-btn:hidden').stop(true, true).fadeIn();
		} else {
			$('#scroll-up-btn').stop(true, true).fadeOut();
		}	
	});
	$('#scroll-up-btn').on('click',function(){
		$('html, body').stop().animate({
        		scrollTop: 0
    		}, 300, 'linear');
		return false;
	});

//////////////////////////////////////
//sample metadata
	checkSampleMetadata();
	
	$( "#human-gender" ).fadeOut('fast');
	$( "#human-age" ).fadeOut('fast');

	$( '#pg-cb-gender' ).on("change",function(){
		if($(this).is(':checked')){
			$( "#human-gender" ).fadeIn('fast');
		} else {
			$( "#human-gender" ).fadeOut('fast');
		}
	});
	$( '#pg-cb-age' ).on("change",function(){
		if($(this).is(':checked')){
			$( "#human-age" ).fadeIn('fast');
		} else {
			$( "#human-age" ).fadeOut('fast');
		}
	});


	//sample type
	if($( '#edge-sample-type').val() == "human") {
			$( "#edge-pg-human" ).fadeIn('fast');
			$( "#edge-pg-host-associated" ).fadeIn('fast');
			$( "#edge-sample-source-options-host" ).fadeIn('fast');
			$( "#edge-sample-source-options-nonhost" ).fadeOut('fast');
	}
	if($( '#edge-sample-type' ).val() == "animal") {
			$( "#edge-pg-human" ).fadeOut('fast');
			$( "#edge-pg-host-associated" ).fadeIn('fast');
			$( "#edge-sample-source-options-host" ).fadeIn('fast');
			$( "#edge-sample-source-options-nonhost" ).fadeOut('fast');
	}
	if($( '#edge-sample-type' ).val() == "environmental") {
			$( "#edge-pg-human" ).fadeOut('fast');
			$( "#edge-pg-host-associated" ).fadeOut('fast');
			$( "#edge-sample-source-options-host" ).fadeOut('fast');
			$( "#edge-sample-source-options-nonhost" ).fadeIn('fast');
	}

	$( '#edge-sample-type' ).on("change",function(){
		if ( $(this).val() == 'human' ){
			$( "#edge-pg-human" ).fadeIn('fast');
			$( "#edge-pg-host-associated" ).fadeIn('fast');
			$('#edge-pg-host').val("human");
			$( "#edge-sample-source-options-host" ).fadeIn('fast');
			$( "#edge-sample-source-options-nonhost" ).fadeOut('fast');

		}
		else if ($(this).val() == 'environmental' ){
			$( "#edge-pg-human" ).fadeOut('fast');
			$( "#edge-pg-host-associated" ).fadeOut('fast');
			$( "#edge-sample-source-options-host" ).fadeOut('fast');
			$( "#edge-sample-source-options-nonhost" ).fadeIn('fast');
			
		}
		else if ($(this).val() == 'animal' ){
			$( "#edge-pg-human" ).fadeOut('fast');
			$( "#edge-pg-host-associated" ).fadeIn('fast');
			$('#edge-pg-host').val("");
			$( "#edge-sample-source-options-host" ).fadeIn('fast');
			$( "#edge-sample-source-options-nonhost" ).fadeOut('fast');
		}
	});

	//seq platform
	if($( '#edge-pg-seq-platform' ).val() =='Illumina') {
			$( "#edge-pg-sequencer-options-ill" ).fadeIn('fast');
			$( "#edge-pg-sequencer-options-ion" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-nan" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-pac" ).fadeOut('fast');

	}
	if($( '#edge-pg-seq-platform' ).val() =='IonTorrent') {
			$( "#edge-pg-sequencer-options-ill" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-ion" ).fadeIn('fast');
			$( "#edge-pg-sequencer-options-nan" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-pac" ).fadeOut('fast');

	}
	if($( '#edge-pg-seq-platform' ).val() =='Nanopore') {
			$( "#edge-pg-sequencer-options-ill" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-ion" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-nan" ).fadeIn('fast');
			$( "#edge-pg-sequencer-options-pac" ).fadeOut('fast');

	}
	if($( '#edge-pg-seq-platform' ).val() =='PacBio') {
			$( "#edge-pg-sequencer-options-ill" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-ion" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-nan" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-pac" ).fadeIn('fast');

	}

	$( '#edge-pg-seq-platform' ).on("change",function(){
		if ( $(this).val() == 'Illumina' ){
			$( "#edge-pg-sequencer-options-ill" ).fadeIn('fast');
			$( "#edge-pg-sequencer-options-ion" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-nan" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-pac" ).fadeOut('fast');
		}
		else if ($(this).val() == 'IonTorrent' ){
			$( "#edge-pg-sequencer-options-ill" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-ion" ).fadeIn('fast');
			$( "#edge-pg-sequencer-options-nan" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-pac" ).fadeOut('fast');
		}
		else if ($(this).val() == 'Nanopore' ){
			$( "#edge-pg-sequencer-options-ill" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-ion" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-nan" ).fadeIn('fast');
			$( "#edge-pg-sequencer-options-pac" ).fadeOut('fast');
		}
		else if ($(this).val() == 'PacBio' ){
			$( "#edge-pg-sequencer-options-ill" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-ion" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-nan" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-pac" ).fadeIn('fast');
		}
	});

	//update sample metadata
	$("#edge-pg-metadata-edit").on("click", function(){
		page.find( ".edge-action-panel" ).panel( "close" );
		editSampleMetadata();
	});



	function checkSampleMetadata(){
    		$.ajax({
    			type: "POST",
			url: "./cgi-bin/edge_sample_metadata.cgi", 
        		cache: false,
        		dataType: "json",
        		data: {"action": 'check' },
				error: function(XMLHttpRequest, textStatus, errorThrown) { 
         			console.log("ERROR");
        		}, // error 
        		// script call was successful 
        		// data contains the JSON values returned by the cgi script 
        		success: function(data){
				console.log("sample metadata is " + data.metadata);
				if (data.metadata == "on") { // sample_metadata_bsve is on
					//console.log("is on");
					$( "#edge-sample-metadata" ).fadeIn('fast');
					localStorage.metadata = true;
          			} else { 
					//console.log("is off");
					$( "#edge-sample-metadata" ).fadeOut('fast');
					localStorage.metadata = false;
				}	
       			} // success
    		});
	}

	function editSampleMetadata() {
		if (!focusProjName){
			showWarning("No focus project for the action");
			return false;
		}
		$.ajax({
			url: "./cgi-bin/edge_sample_metadata.cgi",
			type: "POST",
			dataType: "html",
			cache: false,
			data: { "proj" : focusProjName, "projname":focusProjRealName, "sid":localStorage.sid, "action":"edit" },
			beforeSend: function(){
				$.mobile.loading( "show", {
					text: "Load...",
					textVisible: 1,
					html: ""
				});
			},
			complete: function() {
				$.mobile.loading( "hide" );
			},
			success: function(data){
				//console.log("got response");
				allMainPage.hide();
				$( "#edge-content-report" ).html(data);
				$( "#edge-content-report div[data-role='popup']" ).popup();
				$( "#edge-content-report > div[data-role='collapsible'] table " ).table();
				$( "#edge-content-report > div[data-role='collapsible']" ).collapsible();
				$( "#edge-content-report fieldset[data-role='controlgroup']" ).controlgroup();
				$( "#edge-content-report" ).show();
				$( "#edge-content-report" ).find("img").lazyLoadXT();
				$( "#edge-content-report" ).find("iframe").lazyLoadXT();
				$( "#edge-content-report" ).enhanceWithin();
								
				$.getScript( "./javascript/edge-metadata.js" )
					.done(function( script, textStatus ) {
					//	console.log( "edge-output.js loaded: " + textStatus );
					})
					.fail(function( jqxhr, settings, exception ) {
						console.log( jqxhr, settings, exception );
					});
				
			},
			error: function(data){
				$.mobile.loading( "hide" );
				showWarning("Failed to retrieve the report. Please REFRESH the page and try again."+data);
			}
		});
	};
	
	//init Host list
	$.getJSON( "data/host_list.json", function( data ) {
                var genomeIds = Object.keys(data);
                genomeIds.sort();

                $.each( genomeIds, function( key, val ) {
                        var name=val;
                        name = name.replace(/_/g, " ");
                        $("#edge-hostrm-file-fromlist").append($("<option value="+val+">"+name+"</option>"));
                });
                $("#edge-hostrm-file-fromlist").selectmenu( "refresh" );
        });

	//geo location
        $("#geocomplete").geocomplete({
          map: ".map_canvas",
          details: "form",
          types: ["geocode", "establishment"],
        });

//END sample metatdata	
});
