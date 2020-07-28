$( document ).ready(function()
{
	var page = $( this );
	//with checkbox
	//https://www.gyrocode.com/projects/jquery-datatables-checkboxes/
	var ProjDataTable = $('#edge-gisaid-metadata-project-page-table').DataTable({
			"pageLength": 25,
			"columnDefs": [ 
				{'targets': 0, 'checkboxes': {'selectRow': true} },
				{'targets':[2,3,4,5,6,7,8,9,10], 'type': 'string', 'width': "20%",
					'render': function (data, type, row, meta) {
						if (type === 'filter' || type === 'sort') {
							var api = new $.fn.dataTable.Api(meta.settings);
							var td = api.cell({row: meta.row, column: meta.col}).node();
							data = $('select, input[type="text"]', td).val();
						}
						return data;
					}
				}
 			],
			"select": { "style": 'multi'},
			"scrollX": true,
			"destroy": true,
			"order": [[ 1, 'asc' ]],
			"initComplete": function() {
				$(this.api().table().container()).find('input[type="search"]').parent().wrap('<form>').parent().attr('autocomplete','off').css('overflow','hidden').css('margin','auto');
				if ( $('.dt-checkboxes').length === 1 ){
					this.api().column(0).checkboxes.select();
				}
			},
			"drawCallback" : function(settings){
				$( ".edge-project-page-link").unbind('click').on('click', function(e){
					e.preventDefault();
					var pname = $(this).attr("data-pid");
					updateReport(pname);
				});
			},
			"rowCallback": function( nRow, aData, iDisplayIndex ) {
			}
	});
	// For search/filter input and select datatable
	//https://www.gyrocode.com/articles/jquery-datatables-how-to-search-and-order-by-input-or-select-elements/
	$( '#edge-gisaid-metadata-project-page-table' ).on( 'change', 'tbody select, tbody input[type="text"]', function () {
 		//invalidate the DT cache
		ProjDataTable.cell($(this).closest('td')).invalidate();            
	} );
	// update/download batch form
	$("#edge-gisaid-form-batch-update,#edge-gisaid-form-batch-download").on( "click", function() {
		var rows_selected_projCodes = ProjDataTable.column(0).checkboxes.selected();
		var projs = rows_selected_projCodes.join();
		if (rows_selected_projCodes.length === 0 ){
			showWarning("There are no projects selected.");
			return;
		}
		//console.log(rows_selected_projCodes);
		var action = (this.id.toLowerCase().indexOf("download") >= 0)? "batch-download": (this.id.toLowerCase().indexOf("update") >= 0)? "batch-update": "batch-upload2gisaid";
		var formDom = $("#edge-gisaid-batch-upload-form");
		gisaid_actions(projs,formDom,action);
	});
	( navigator.appVersion.indexOf("Mac")>=0)? $("#mac-scroll-bar-note").show():$("#mac-scroll-bar-note").hide();
	
	// cancel batch form
	$("#edge-gisaid-form-batch-cancel").on( "click", function() {
		$( "#edge-gisaid-metadata-project-page" ).hide();
		$( "#edge-project-page" ).show();
	});

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
		var formDom = $("#edge-gisaid-upload-form");
		gisaid_actions(proj,formDom,action);
		
	});

	function gisaid_actions(proj, form, action){
		$.ajax({
				url: "./cgi-bin/edge_gisaid_upload.cgi",
				type: "POST",
				dataType: "json",
				cache: false,
				//data: $( "#edge-run-pipeline-form" ).serialize(),
				data: ( form.serialize() +'&'+ $.param({ 'proj': proj, 'sid':localStorage.sid, 'action':action, "userDir":localStorage.udir })),
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
						if ( action.toLowerCase().indexOf("download") >= 0){
							window.open(edge_path + obj.PATH);
						}
						if ( action.toLowerCase().indexOf("upload2gisaid") >= 0){
							$( "#edge-submit-info" ).fadeIn("fast");
							//show message
							$( "#edge_integrity_dialog_header" ).text("Message");
							$( "#edge_integrity_dialog_content" ).text("Your project was successfully submitted to the GISAID.");
							setTimeout( function() { $( "#edge_integrity_dialog" ).popup('open'); }, 300 );
							if (action == "batch-upload2gisaid"){
								$( "#edge-gisaid-metadata-project-page" ).hide();
								$( "#edge-project-page" ).show();
							}else{
								updateReport($('#edge-output-projid').attr("data-pid"));
							}
						}
						if (action.toLowerCase().indexOf("update") >= 0){
							$( "#edge_integrity_dialog_header" ).text("Message");
							$( "#edge_integrity_dialog_content" ).text("Your project(s) metadata was successfully updated.");
							setTimeout( function() { $( "#edge_integrity_dialog" ).popup('open'); }, 300 );
							if (action == "batch-update"){
								$( "#edge-gisaid-metadata-project-page" ).hide();
								$( "#edge-project-page" ).show();
							}else{
								updateReport($('#edge-output-projid').attr("data-pid"));
							}
						}
					}
					else{
						// display error information
						if(! $.isEmptyObject(obj.PARAMS)){
							var projlistIndex = Object.keys(obj.PARAMS);
							$.each(projlistIndex, function(j,h){
								$.each(obj.PARAMS[h], function(i,v){
									$("[name^='"+i+"']").eq(h).addClass("highlight");
									$("#"+i).addClass("highlight");
									$("#"+i).parents('div[data-role="collapsible"]').collapsible( "option", "collapsed", false );
									$( "#edge-submit-info" ).fadeIn("fast");
									var dom = "<li data-icon='delete' data-theme='c' class='list-info-delete'><a href='#'>"+v+"</a></li>";
							
									$(dom).appendTo("#edge-submit-info").on("click", function(){
										$('html, body').animate({
											scrollTop: $("#"+i).offset().top-100
										}, 200);
									});
								});
							});
						}else{
							$( "#edge-submit-info" ).fadeIn("fast");
							var dom;
							if ( action.toLowerCase().indexOf("upload2gisaid") >= 0){
								dom = "<li data-icon='delete' data-theme='c' class='list-info-delete'><a href='#'>FAILED to submit to GISAID. Please check GISAID/submit.log in the project directory for detail.</a></li>";
							}
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
	}
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
