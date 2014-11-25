$( document ).ready(function()
{
	var running_proj=0;
	var finished_proj=0;
	var failed_proj=0;
	var focusProjName;
	var focusProjStatus;
	var focusProjTime;
	var focusProjRunningStep;

	var interval = 1000*5; //check every 7 secs
	var updateProjInterval;
	var updateLogInterval;
	var inputFileID;
	var inputFileDir  = "/";
	
	var page = $( this );
	var allMainPage = $(".edge-main-page");

	//init page
	sync_input(); //sync input with switch
	updateProject();
	integrityCheck();
	
	allMainPage.hide();
	$( "#edge-content-home" ).fadeIn();
	
	$( "a[href=#edge-content-home]" ).on( "click", function(){
		allMainPage.hide();
		$( "#edge-content-home" ).fadeIn();
		page.find( ".edge-navmenu-panel:not(.edge-panel-page-nav)" ).panel( "close" );
	});
	$( "a[href=#edge-content-pipeline]" ).on( "click", function(){
		allMainPage.hide();
		$( "#edge-submit-info" ).children().remove();
		$( "#edge-content-pipeline" ).fadeIn();
		page.find( ".edge-navmenu-panel:not(.edge-panel-page-nav)" ).panel( "close" );
	});

	//input sequence
	$( "#edge-input-options" ).hide();
	$( "#edge-input-toggle" ).on( "click", function(){
		$( "#edge-input-options" ).toggle();
	});
	
    	// batch input 
    	$('#edge-batch-sample-input').click( function(e) {
    		e.preventDefault();
    		var sampleInput = "#each unique project name in the bracket []\n" + "[Project1]\n" + "#q1=/path/to/paired_end_file_1\n" + "q1=edgeui_input/testData/Ecoli_10x.1.fastq\n" + 
    	                  "#q2=/path/to/paired_end_file_2\n" + "q2=edgeui_input/testData/Ecoli_10x.2.fastq\n" + "description=\"test batch input project 1\"\n";
    		sampleInput = sampleInput + "[Project2]\n" + "q1=edgeui_input/testData/Ecoli_10x.1.fastq\n" + "q2=edgeui_input/testData/Ecoli_10x.2.fastq\n" + "description=\"test batch input project 2\"\n";
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
			$( "#edge_integrity_dialog_content" ).text(f.name + " is not a valid text file.");
			$( "#edge_integrity_dialog" ).popup('open');
		} else if (f.size > 5120) {
	         	// > 5Kb  It seems too big to process on the text 
			$( "#edge_integrity_dialog_content" ).text(f.name + " with size " + f.size + " bytes is too big (>5k)");
			$( "#edge_integrity_dialog" ).popup('open');
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
  	
	//init reserch usage bar
	$( "input[id$='usage-bar']" ).each(function( event ) {
		$(this).parent().find("input").hide();
		$(this).parent().find(".ui-slider-handle").remove();
		$(this).parent().find(".ui-slider-track").css('margin','0 15px 0 15px').css('pointer-events','none');
	});

	//init file tree
	$( "#edge_file_tree" ).fileTree({
			root: inputFileDir,
			script: './cgi-bin/jqueryFileTree.cgi',
		}, function(file) {
			$( "#"+inputFileID ).val(file);
			$( "#edge_file_dialog" ).popup('close');
		}
	);
 
	//show file browser
	$( ".edge-file-selector" ).on( "click", function() {
		inputFileID = $(this).prevAll().children().prop("id");
	});
	
	// bind event to content
	$( ".ui-btn.ui-input-clear" ).on( "click", function() {
		var dvalue = $(this).prev("input").prop( 'defaultValue' );
		$(this).prev("input").val(dvalue);
	});

	//button for adding input fileds
    	$('#btnAdd-edge-input-se').click( function(e) {
		e.preventDefault();
        	// how many "duplicatable" input fields we currently have
        	var num = $('.edge-input-se-block').length;	
        	
        	// the numeric ID of the new input field being added	
        	var newNum	= new Number(num + 1);		
        	var newElem = $('#edge-input-se-block' + num ).clone().attr('id', 'edge-input-se-block' + newNum);
        	newElem.find('label').attr( 'for', 'edge-input-se' + newNum ).text('Single-end FASTQ file (' + newNum + ')');            
        	newElem.find('input').attr( 'id', 'edge-input-se' + newNum ).attr('name', 'edge-input-se[]');
        	newElem.find('.btnDel-edge-input-se').css("visibility","visible");
        	// insert newElem
        	$('#edge-input-se-block' + num).after(newElem);
        	
        	// bind the selector 
        	newElem.find(".edge-file-selector").on( "click", function() {
				inputFileID = 'edge-input-se' + newNum;
		    });
		    
		    newElem.find(".btnDel-edge-input-se").on( "click", function() {
				$('#edge-input-se-block' + newNum ).remove();
				$('#btnAdd-edge-input-se' ).removeClass('ui-disabled');
		    });
        	// business rule: limit the number of fields to 5
        	if (newNum == 5) {
        	     $('#btnAdd-edge-input-se' ).addClass('ui-disabled');
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
  
	//init host
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

	$.mobile.document
    .on( "listviewcreate", "#filter-menu-menu,#edge-hostrm-file-fromlist-menu", function( event ) {
        var input,
            list = $( event.target ),
            form = list.jqmData( "filter-form" );
        if ( !form ) {
            input = $( "<input data-type='search'></input>" );
            form = $( "<form></form>" ).append( input );
            input.textinput();
            list
                .before( form )
                .jqmData( "filter-form", form ) ;
            form.jqmData( "listview", list );
        }
        list.filterable({
            input: input,
            children: "> li:not(:jqmData(placeholder='true'))"
        });
    })
    .on( "pagecontainerbeforeshow", function( event, data ) {
        var listview, form,
            id = data.toPage && data.toPage.attr( "id" );
        if ( !( id === "filter-menu-dialog" || id === "edge-hostrm-file-fromlist-dialog" ) ) {
            return;
        }
        listview = data.toPage.find( "ul" );
        form = listview.jqmData( "filter-form" );
        data.toPage.jqmData( "listview", listview );
        listview.before( form );
    })
    .on( "pagecontainerhide", function( event, data ) {
        var listview, form,
            id = data.toPage && data.toPage.attr( "id" );
        if ( !( id === "filter-menu-dialog" || id === "edge-hostrm-file-fromlist-dialog" ) ) {
            return;
        }
        listview = data.toPage.jqmData( "listview" ),
        form = listview.jqmData( "filter-form" );
        listview.before( form );
    });


	//init ref
	$.getJSON( "data/ref_list.json", function( data ) {
		var genomeIds = Object.keys(data);
		genomeIds.sort();

		$.each( genomeIds, function( key, val ) {
			var name=val;
			name = name.replace(/_uid\d+$/, "");
			name = name.replace(/_/g, " ");
			$("#edge-ref-file-fromlist").append($("<option value="+val+">"+name+"</option>"));
		});
		$("#edge-ref-file-fromlist").selectmenu( "refresh" );
	});

	$.mobile.document
    .on( "listviewcreate", "#filter-menu-menu,#edge-ref-file-fromlist-menu", function( event ) {
        var input,
            list = $( event.target ),
            form = list.jqmData( "filter-form" );
        if ( !form ) {
            input = $( "<input data-type='search'></input>" );
            form = $( "<form></form>" ).append( input );
            input.textinput();
            list
                .before( form )
                .jqmData( "filter-form", form ) ;
            form.jqmData( "listview", list );
        }
        list.filterable({
            input: input,
            children: "> li:not(:jqmData(placeholder='true'))"
        });
    })
    .on( "pagecontainerbeforeshow", function( event, data ) {
        var listview, form,
            id = data.toPage && data.toPage.attr( "id" );
        if ( !( id === "filter-menu-dialog" || id === "edge-ref-file-fromlist-dialog" ) ) {
            return;
        }
        listview = data.toPage.find( "ul" );
        form = listview.jqmData( "filter-form" );
        data.toPage.jqmData( "listview", listview );
        listview.before( form );
    })
    .on( "pagecontainerhide", function( event, data ) {
        var listview, form,
            id = data.toPage && data.toPage.attr( "id" );
        if ( !( id === "filter-menu-dialog" || id === "edge-ref-file-fromlist-dialog" ) ) {
            return;
        }
        listview = data.toPage.jqmData( "listview" ),
        form = listview.jqmData( "filter-form" );
        listview.before( form );
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
		sync_input();
	});

	$( ".edge-collapsible-options div.ui-slider-switch" ).on( "mouseover", function() {
		$(this).parents('div[data-role="collapsible"]').collapsible("disable");
	});
	$( ".edge-collapsible-options div.ui-slider-switch" ).on( "mouseout", function() {
		$(this).parents('div[data-role="collapsible"]').collapsible("enable");
	});

	$( "#edge-view-log-btn" ).on( "mouseover", function() {
		getLog("./EDGE_output/"+focusProjName+"/process_current.log");
	});

	$( "#edge_log_dialog" ).popup({
		beforeposition: function( event, ui ) {
			getLog("./EDGE_output/"+focusProjName+"/process_current.log");
		},
		afteropen: function( event, ui ) {
			if( focusProjStatus != "finished" ){
				updateLogInterval = setInterval(
					function(){
						getLog("./EDGE_output/"+focusProjName+"/process_current.log");
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

	//confirm dialog
	$( "a[id^='action']" ).on('mouseover', function(){
		var action = $(this).attr("data");
		$("#edge_confirm_dialog_content").html("Do you want to <span id='action_type'>"+action.toUpperCase()+"</span> project "+focusProjName+"?");
	});


	$( "#edge_confirm_dialog" ).popup({
		afterclose: function( event, ui ) {
			updateProject(focusProjName);
		}
	});

	$("#edge_confirm_dialog a:contains('Confirm')").on("click",function(){
		var action = $("#edge_confirm_dialog_content > span").html();
		$.ajax({
			url: "./cgi-bin/edge_action.cgi",
			type: "POST",
			dataType: "json",
			cache: false,
			data: { "proj" : focusProjName, "action": action },
			beforeSend: function(){
				$.mobile.loading( "show", {
					text: "Executing "+ action.toUpperCase() +" command...",
					textVisible: 1,
					html: ""
				});
			},
			complete: function() {
				$.mobile.loading( "hide" );
				updateProject(focusProjName);
			},
			success: function(data){
				if( data.STATUS == "SUCCESS" ){
					$.mobile.loading( "hide" );
					$( "#edge_integrity_dialog_content" ).text( data.INFO );
					$( "#edge_integrity_dialog" ).popup('open');
				}
				else{
					$.mobile.loading( "hide" );
					$( "#edge_integrity_dialog_content" ).text(data.INFO);
					$( "#edge_integrity_dialog" ).popup('open');
				}
			},
			error: function(data){
				$.mobile.loading( "hide" );
				$( "#edge_integrity_dialog_content" ).text("ACTION FAILED: Please try again or contact your system administrator.");
				$( "#edge_integrity_dialog" ).popup('open');
			}
		});
	});

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

	//update report
	function updateReport(pname) {
		$.ajax({
			url: "./cgi-bin/edge_report.cgi",
			type: "POST",
			dataType: "html",
			cache: false,
			data: { "proj" : pname },
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
				allMainPage.hide();
				$( "#edge-content-report" ).html(data);
				$( "#edge-content-report div[data-role='popup']" ).popup();
				$( "#edge-content-report > div[data-role='collapsible'] table " ).table();
				$( "#edge-content-report > div[data-role='collapsible']" ).collapsible();
				$( "#edge-content-report" ).show();

				$( "#edge-content-report" ).find("img").lazyLoadXT();
				$( "#edge-content-report" ).find("iframe").lazyLoadXT();
				
				$.getScript( "./javascript/edge-output.js" )
					.done(function( script, textStatus ) {
						console.log( "edge-output.js loaded: ".textStatus );
					})
					.fail(function( jqxhr, settings, exception ) {
						console.log( "edge-output.js loaded: ".exception );
					});

				var pstatus = $( "#edge-content-report" ).find("p:first" ).text().match( /Project Status: (.+)/m );

				if( pstatus[1] != "Complete" ){
					$( "#edge-content-report div.ui-grid-a" ).hide();
					$( "#edge-content-report div.ui-grid-c" ).hide();
				}
			},
			error: function(data){
				$.mobile.loading( "hide" );
				$( "#edge_integrity_dialog_content" ).text("Failed to retrieve the report. Please REFRESH the page and try again.");
				$( "#edge_integrity_dialog" ).popup('open');
			}
		});
	}

	//submit
	$( "#edge-form-submit" ).on( "click", function() {
		$.ajax({
			url: "./cgi-bin/edge_submit.cgi",
			type: "POST",
			dataType: "json",
			cache: false,
			data: page.find("form").serialize(),
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
				updateProject();
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
						}
						$( "#edge-submit-info" ).append(dom).fadeIn("fast");
					}
				});

				if( obj.SUBMISSION_STATUS == "success" ){
					$( "#edge-submit-info" ).fadeIn("fast");
					var dom = "<li data-icon='info' class='list-info'><a href='#'>The job has been submitted successfully. Click to see open pregress panel.</a></li>";
					
					$(dom).on( "click", function(){page.find( ".edge-action-panel" ).panel( "open" );})
					      .appendTo( "#edge-submit-info");
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
		$("#edge-anno-sw").on("change", function(){
			if ( $(this).val() == 1 ){
				if ( $( "#edge-assembly-sw" ).val() == 0 ){
					$( "#edge_integrity_dialog_content" ).text("Annotation function can only be performed when assembly function turns on.");
					$( "#edge_integrity_dialog" ).popup('open');
					$( "#edge-anno-sw" ).val(0).slider("refresh");
				}
			}
		});
		$("#edge-primer-valid-sw").on("change", function(){
			if ( $(this).val() == 1 ){
				if ( $( "#edge-assembly-sw" ).val() == 0 ){
					$( "#edge_integrity_dialog_content" ).text("Primer validation function can only be performed when assembly function turns on.");
					$( "#edge_integrity_dialog" ).popup('open');
					$( "#edge-primer-valid-sw" ).val(0).slider("refresh");
				}
			}
		});
		$("#edge-primer-adj-sw").on("change", function(){
			if ( $(this).val() == 1 ){
				if ( $( "#edge-assembly-sw" ).val() == 0 ){
					$( "#edge_integrity_dialog_content" ).text("Primer design function can only be performed when assembly function turns on.");
					$( "#edge_integrity_dialog" ).popup('open');
					$( "#edge-primer-adj-sw" ).val(0).slider("refresh");
				}
			}
		});

		$("#edge-assembly-sw").on("change", function(){
			if ( $(this).val() == 0 ){
				if ( $( "#edge-anno-sw" ).val() == 1 ){
					$( "#edge_integrity_dialog_content" ).text("Assembly function can not be turned off because annotation function is on.");
					$( "#edge_integrity_dialog" ).popup('open');
					$( "#edge-assembly-sw" ).val(1).slider("refresh");
				}
				else if ( $( "#edge-primer-valid-sw" ).val() == 1 ){
					$( "#edge_integrity_dialog_content" ).text("Assembly function can not be turned off because primer validation function is on.");
					$( "#edge_integrity_dialog" ).popup('open');
					$( "#edge-assembly-sw" ).val(1).slider("refresh");
				}
				else if ( $( "#edge-primer-adj-sw" ).val() == 1 ){
					$( "#edge_integrity_dialog_content" ).text("Assembly function can not be turned off because primer design function is on.");
					$( "#edge_integrity_dialog" ).popup('open');
					$( "#edge-assembly-sw" ).val(1).slider("refresh");
				}
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

	function updateProject(pname) {
		focusProjName = pname;
		$.ajax({
			url: './cgi-bin/edge_info.cgi',
			type: "POST",
			dataType: "json",
			cache: false,
			data: { "proj" : pname },
			complete: function(data){
				console.log("finished_proj="+finished_proj);
				console.log("running_proj="+running_proj);
				console.log("failed_proj="+failed_proj);
				console.log("focusProjName="+focusProjName);
				console.log("focusProjStatus="+focusProjStatus);
				$( "#edge-project-list-ul > li" ).on("click", function(){
					var pname = $(this).children("a").text();
					updateProject(pname);
					updateReport(pname);
				});

				clearInterval(updateProjInterval);
				updateProjInterval = setInterval(function(){updateProject(pname)}, interval);
			},
			success: function(obj) {
				// update focused proj
				focusProjName   = obj.INFO.NAME;
				focusProjStatus = obj.INFO.STATUS;
				focusProjTime   = obj.INFO.TIME;
				finished_proj   = 0;
				running_proj    = 0;
				failed_proj     = 0;

				//resource usage bar
				$("#cpu-usage-bar").val(obj.INFO.CPUU).slider("refresh");
				$("#mem-usage-bar").val(obj.INFO.MEMU).slider("refresh");
				$("#disk-usage-bar").val(obj.INFO.DISKU).slider("refresh");
				$("#cpu-usage-val").html(obj.INFO.CPUU+" %");
				$("#mem-usage-val").html(obj.INFO.MEMU+" %");
				$("#disk-usage-val").html(obj.INFO.DISKU+" %");

				// project list
				if(! $.isEmptyObject(obj.LIST)){
					$( "#edge-project-list-ul .edge-proj-list-li" ).remove();
		
					var listIdOrder = Object.keys(obj.LIST);
					listIdOrder.sort(function(a,b){ return obj.LIST[a].TIME < obj.LIST[b].TIME ? -1 : obj.LIST[a].TIME > obj.LIST[b].TIME; }).reverse();

					$.each(listIdOrder, function(i,v){
						var proj_list_obj = obj.LIST[v];
						if( proj_list_obj.NAME ){
							var name   = proj_list_obj.NAME;
							var time   = proj_list_obj.TIME;
							var pstatus = proj_list_obj.STATUS;
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
  								case "failed":
									projClass = "edge-time-bg-red";
									projIcon  = "ui-icon-delete";
									failed_proj++;
									break;
  								default:
  									projClass = "edge-time-bg-grey";
  									projIcon  = "ui-icon-refresh";
							}

							var dom = "<li class='edge-proj-list-li'><div class='edge-project-time "+projClass+"'>"+time+"</div><a href='' class='edge-project-list ui-btn ui-btn-icon-right "+projIcon+"' title='"+desc+"'>"+name+"</a></li>";
							$(dom).appendTo( "#edge-project-list-ul" );
						}
					});
				}
				
				if( $( ".edge-proj-list-li" ).size() == 0 ){
					var dom = "<li class='edge-proj-list-li'><a href='#' class='edge-project-list ui-btn ui-btn-icon-right ui-icon-check'>No project found</a></li>";
					$( "#edge-project-list-ul" ).append(dom);
				}

				// progress info
				if(! $.isEmptyObject(obj.PROG))
				{
					$( "#edge-progress-ul > li" ).fadeOut().remove();

					var dom = "<li data-role='list-divider'>"+obj.INFO.NAME+"</li>";
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
								return true;
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
						$( "#edge-progress-ul" ).append(dom);
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
				}
			},
			error: function(data){
				$( "#edge-submit-info" ).fadeIn("fast");
				var dom = "<li data-icon='delete' data-theme='c' class='list-info-delete'><a href='#'>FAILED to retrieve project info. Please check server error log for detail.</a></li>";
				$( "#edge-submit-info" ).append(dom);
			}
		});
	};
});
