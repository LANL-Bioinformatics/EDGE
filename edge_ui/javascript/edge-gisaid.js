$( document ).ready(function()
{
	var page = $( this );
	//datepicker
	$('.metadata-input-date').datepicker({
		changeMonth: true,
		changeYear: true,
		dateFormat: 'yy-mm-dd'
	});

	//form reset
	$( "#edge-gisaid-form-reset" ).on( "click", function() {
		$( "#edge-gisaid-upload-form" )[0].reset();
	});

	//form cancel
	$( "#edge-gisaid-form-cancel" ).on( "click", function() {
		updateReport($('#edge-output-projid').attr("data-pid"));
	});

	var loc = window.location.pathname.replace("//","/");
        var edge_path = loc.substring(0,loc.lastIndexOf('/'));
	//$("#edge-gisaid-form-submit").parent().hide();
	//$("#edge-gisaid-form-download").parent().hide();
	//form submit
	//$( "#edge-gisaid-form-submit" ).on( "click", function() {
	$( "#edge-gisaid-form-submit,#edge-gisaid-form-download, #edge-gisaid-form-update" ).on( "click", function() {
		var action = (this.id.toLowerCase().indexOf("download") >= 0)? "download": (this.id.toLowerCase().indexOf("update") >= 0)? "update": "upload2gisaid";
		var proj=$('#edge-output-projid').attr("data-pid");
		$.ajax({
			url: "./cgi-bin/edge_gisaid_upload.cgi",
			type: "POST",
			dataType: "json",
			cache: false,
			//data: $( "#edge-run-pipeline-form" ).serialize(),
			data: ( $("#edge-gisaid-upload-form").serialize() +'&'+ $.param({ 'proj': proj, 'sid':localStorage.sid, 'action':action, "userDir":localStorage.udir })),
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
					if( i!="PARAMS" && i!="SUBMISSION_STATUS" && i!="PATH"){
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
					if ( action == "download"){
						window.open(edge_path + obj.PATH);
					}
					if ( action == "upload2gisaid"){
						$( "#edge-submit-info" ).fadeIn("fast");
						//show message
						$( "#edge_integrity_dialog_header" ).text("Message");
						$( "#edge_integrity_dialog_content" ).text("Your project was successfully submitted to the GISAID.");
						setTimeout( function() { $( "#edge_integrity_dialog" ).popup('open'); }, 300 );

						updateReport($('#edge-output-projid').attr("data-pid"));
					}
					if (action == "update"){
						$( "#edge_integrity_dialog_header" ).text("Message");
						$( "#edge_integrity_dialog_content" ).text("Your project metadata was successfully updated.");
						setTimeout( function() { $( "#edge_integrity_dialog" ).popup('open'); }, 300 );
						updateReport($('#edge-output-projid').attr("data-pid"));
					}
				}
				else{
					// display error information
					if(! $.isEmptyObject(obj.PARAMS)){
						$.each(obj.PARAMS, function(i,v){
							$("#"+i).addClass("highlight");
							$("#"+i).parents('div[data-role="collapsible"]').collapsible( "option", "collapsed", false );
							$( "#edge-submit-info" ).fadeIn("fast");
							var label = $("label[for='"+i+"']").text();
							var dom = "<li data-icon='delete' data-theme='c' class='list-info-delete'><a href='#'>"+v+"</a></li>";
							
							$(dom).appendTo("#edge-submit-info").on("click", function(){
								$('html, body').animate({
									scrollTop: $("#"+i).offset().top-80
								}, 200);
							});
						});
					}else{
						$( "#edge-submit-info" ).fadeIn("fast");
						var dom = "<li data-icon='delete' data-theme='c' class='list-info-delete'><a href='#'>FAILED to submit to GISAID. Please check GISAID/submit.log in the project directory for detail.</a></li>";
						$( "#edge-submit-info" ).append(dom).fadeIn("fast");
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
				var dom = "<li data-icon='delete' data-theme='c' class='list-info-delete'><a href='#'>FAILED to run gisaid submission CGI. Please check server error log for detail.</a></li>";
				$( "#edge-submit-info" ).append(dom).fadeIn("fast");
			}
		});
	});

	function updateReport(pname) {
		$.ajax({
			url: "./cgi-bin/edge_report.cgi",
			type: "POST",
			dataType: "html",
			cache: false,
			data: { "proj" : pname, "sid":localStorage.sid },
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
				
				$( "#edge-content-report" ).html(data);
				$( "#edge-content-report div[data-role='popup']" ).popup();
				$( "#edge-content-report > div[data-role='collapsible'] table " ).table();
				$( "#edge-content-report > div[data-role='collapsible']" ).collapsible();
				$( "#edge-content-report fieldset[data-role='controlgroup']" ).controlgroup();
				$( "#edge-content-report" ).show();
				$( "#edge-content-report" ).find("img").lazyLoadXT();
				$( "#edge-content-report" ).find("iframe").lazyLoadXT();
				$( "#edge-content-report" ).enhanceWithin();
								
				$.getScript( "./javascript/edge-output.js" )
					.done(function( script, textStatus ) {
					//	console.log( "edge-output.js loaded: " + textStatus );
					})
					.fail(function( jqxhr, settings, exception ) {
						console.log( jqxhr, settings, exception );
					});
				var pstatus = $( "#edge-content-report" ).find("p:first").text().match( /Project Status: (.+)/m );

				if( pstatus[1] != "Complete" ){
					$( "#edge-content-report div.ui-grid-a" ).hide();
					$( "#edge-content-report div.ui-grid-c" ).hide();
				}
				//console.log($( "#edge-get-contigs-by-taxa" ));
			},
			error: function(data){
				$.mobile.loading( "hide" );
				showWarning("Failed to retrieve the report. Please REFRESH the page and try again."+pname);
			}
		});
	}


	function showWarning( dialog_content ) {
		$( "#edge_integrity_dialog_content" ).html(dialog_content);
		$( "#edge_integrity_dialog" ).popup('open');
	}
});
//# sourceURL=edge-output.js
