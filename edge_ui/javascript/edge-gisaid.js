$( document ).ready(function()
{
	var checkpidInterval;
	var newWindowHeader = "<html><head><title>EDGE COVID-19</title><link rel='stylesheet' href='css/edge-output.css'/></head><div style='background:#02497e;'><h2 style='color:#fff;position:inherit; padding-left:20px;'>EDGE COVID-19</h2></div>";
	var newWindowFooter = "<div id='newWindowSpinner'  class='edge-sp edge-sp-circle'></div><pre style='background-color: rgba(0,0,0,.8);overflow-y: auto;padding: 1em;'><code id='newWindowMsg' style='white-space: pre-wrap;color: white; font: 300 0.8em 'Bitstream Vera Sans Mono', 'Courier', monospace;'></code></pre></body></html>";
        var newTitle = "EDGE COVID-19";
	var page = $( this );
	$("table td, table th, .tooltip").tooltipster({
						multiple:true, maxWidth: '480', interactive: true
					});
	$("#metadata-tab-example").tooltipster(
            'content', $('<span><table border="1" style="font-size:0.8em"><tr><th>project-name</th><th>virus-name</th><th>virus-passage</th><th>sample-collection-date</th><th>sample-location</th><th>sample-host</th><th>sample-gender</th><th>sample-age</th><th>sample-status</th><th>sample-sequencing-tech</th></tr><tr><td>project1</td><td>hCoV-19/USA/NM-LANL-00001/2020</td><td>Original</td><td>2020-09-08</td><td>North America/USA/New Mexico</td><td>Human</td><td>Male</td><td>65</td><td>Live</td><td>illumina</td></tr><tr><td>project2</td><td>hCoV-19/USA/NM-LANL-00002/2020</td><td>Vero</td><td>2020-07-20</td><td>North America/USA/Arizona</td><td>Human</td><td>Female</td><td>50</td><td>unknown</td><td>Nanopore</td></tr></table></span>')
    ).tooltipster('option','maxWidth','700');

	//with checkbox
	//https://www.gyrocode.com/projects/jquery-datatables-checkboxes/
	var ProjDataTable = $('#edge-gisaid-metadata-project-page-table').DataTable({
			"pageLength": 25,
			"columnDefs": [ 
				{'targets': 0, 'checkboxes': {'selectRow': true} },
				{'targets':[2,3,4,5,6,7,8,9,10,11], 'type': 'string', 'width': "20%",
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
				$( '#edge-gisaid-metadata-project-page-table' ).find('select').selectmenu();
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
	// adjust the columns misalignment sometimes.
	ProjDataTable.columns.adjust();
	// For search/filter input and select datatable
	//https://www.gyrocode.com/articles/jquery-datatables-how-to-search-and-order-by-input-or-select-elements/
	$( '#edge-gisaid-metadata-project-page-table' ).on( 'change', 'tbody select, tbody input[type="text"]', function () {
		var $td = $(this).closest('td');
		// update JQM select span selected name
		if ($(this).is("select")){
			$(this).selectmenu('refresh');
	
			//var value = $(this).val() || $td.find('option').eq(0).html();
			//$td.find('span').html(value);
		}
 		//invalidate the DT cache
		ProjDataTable.cell($td).invalidate();

	} );
	// update/download batch form
	$("#edge-gisaid-form-batch-update,#edge-gisaid-form-batch-download,#edge-gisaid-form-batch-submit").on( "click", function() {
		var action="";
		if (this.id.toLowerCase().indexOf("batch-download") >= 0){
			action = "batch-download";
		}else if (this.id.toLowerCase().indexOf("batch-update") >= 0){
			action = "batch-update";
		}else if (this.id.toLowerCase().indexOf("batch-submit") >= 0){ 
			action = "batch-upload2gisaid";
		}
		var formDom = $("#edge-gisaid-batch-upload-form");
		var rows_selected = ProjDataTable.column(0).checkboxes.selected();
		//console.log(rows_selected);
		var projCodes=[];
		var projNames=[];
		var NotReadyCon=[];
		if (rows_selected.length === 0){
			showWarning("There are no projects selected.");
			return;
		}
		rows_selected.each(function (item, index) {
			var rowIdx = item;
			var pname = $(ProjDataTable.cell(rowIdx,1).data()).eq(0).val();
			var pcode = $(ProjDataTable.cell(rowIdx,1).data()).eq(1).val();
			projCodes.push(pcode);
			projNames.push(pname);
			var selectedCon = $(ProjDataTable.cell(rowIdx,11).data()).eq(0).children("option:selected").html() || 
                                          $(ProjDataTable.cell(rowIdx,11).data()).find('span').html() ;
			if (action == 'batch-upload2gisaid' && ! /Ready to Submit/.test(selectedCon)){
				NotReadyCon.push(pname);
			}

		});
		var projs = projCodes.join();
		if (NotReadyCon.length){
			var actionContent = "The consensus genoems of projects:<br/>" + NotReadyCon.join("<br/>") + ' <br/>&nbsp;&nbsp;are NOT Ready to Submit. Do you want to proceed? <p>This action can not be undone.</p>';
			$("#edge_confirm_dialog_content").html(actionContent);
			$('#edge_confirm_dialog').enhanceWithin().popup('open');
			$("#edge_confirm_dialog a:contains('Confirm')").unbind('click').on("click",function(){
				gisaid_actions(projs,formDom,action,projNames.join());
                	});
		}else{
			gisaid_actions(projs,formDom,action,projNames.join());
		}
	});
	$("#edge-gisaid-form-batch-template-download").on( "click", function() {
		action='batch-template-download';
		ProjDataTable.column(0).checkboxes.select();
		var formDom = $("#edge-gisaid-batch-upload-form");
		var rows_selected = ProjDataTable.column(0).checkboxes.selected();
		var projCodes=[];
		var projNames=[];
		rows_selected.each(function (item, index) {
			var rowIdx = item;
			var pcode = $(ProjDataTable.cell(rowIdx,1).data()).eq(1).val();
			var pname = $(ProjDataTable.cell(rowIdx,1).data()).eq(0).val();
			projCodes.push(pcode);
			projNames.push(pname);
		});
		var projs = projCodes.join();
		gisaid_actions(projs,formDom,action,projNames.join());
		ProjDataTable.column(0).checkboxes.deselect();
	});
	
	$("#metadata-upload-file").on('change',function(){
		if (this.files.length){
			Papa.parse(this.files[0], {
				header: true,
				skipEmptyLines: true,
				complete: function(results) {
					// each row 
					$.each(results.data,function(index,item){
						$.each(item, function(key,val){
							if ( $('input[name$="' + key + 's"]').length ){
								$('input[name$="' + key + 's"]').eq(index).val(val);
							}
							if( $('select[name$="' + key + 's"]').length || $('select[name$="' + key + '"]').length  ){
								var res = val.toLowerCase();
								if (!val){
									res = 'unknown';
								}
								if (res != 'unknown'){
									// CAP on first letter
									res = val.charAt(0).toUpperCase() + val.slice(1);
								}
								$('select[name$="' + key + 's"]').eq(index).val(res).selectmenu('refresh');
								$('select[name$="' + key + '"]').eq(index).val(res).selectmenu('refresh');
							} 
						}); 
					});
				}
			});
		}
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
	$( "#edge-sra-form-reset" ).on( "click", function() {
		$( "#edge-sra-upload-form" )[0].reset();
	});
	
	//form cancel
	$( "#edge-gisaid-form-cancel, #edge-sra-form-cancel" ).on( "click", function() {
		updateReport($('#edge-output-projid').attr("data-pid"));
	});
	//template
	$( "#metadata-upload-content" ).hide();
	$( '#metadata-upload-toggle' ).on("click", function() {
		$("#metadata-upload-content").toggle();
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
		var projname=$('#edge-project-title').text().replace(' /','');
		var formDom = $("#edge-gisaid-upload-form");
		var selectedCon = $('#metadata-sample-consensus').children("option:selected").html();
		if ( action == 'upload2gisaid' && ! /Ready to Submit/.test(selectedCon) ){
			var actionContent = selectedCon + ' is NOT Ready to Submit. Do you want to proceed? <p>This action can not be undone.</p>';
			$("#edge_confirm_dialog_content").html(actionContent);
			$('#edge_confirm_dialog').enhanceWithin().popup('open');
			$("#edge_confirm_dialog a:contains('Confirm')").unbind('click').on("click",function(){
				gisaid_actions(proj,formDom,action,projname);
                	});
		}else{
			gisaid_actions(proj,formDom,action,projname);
		}
		
	});

	function gisaid_actions(proj, form, action, projname){
		var w;
		var info_dom_id = ( action.indexOf("batch") >=0 )? "edge-repo-batch-submit-info" : "edge-repo-submit-info";
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
					if (action.toLowerCase().indexOf("upload2gisaid") >= 0){
						w = window.open("","new","width=480,height=480");
						w.document.body.innerHTML = '';
						w.document.write( newWindowHeader + "<p id='newWindowInfo'>Submit project(s) " + projname + " to GISAID and NCBI. Please wait...</p>" + newWindowFooter);
					}
				},
				complete: function(data) {
					$("#" + info_dom_id).listview("refresh");
					setTimeout(function(){ $.mobile.loading( "hide" );},1000);
				},
				success: function(obj){
					// display general submission error
					$.each(obj, function(i,v){
						if( i!="PARAMS" && i!="SUBMISSION_STATUS" && i!="PATH" && i!="PID"){
							var dom;
							if( this.STATUS == "failure" ){
								dom = "<li data-icon='delete' data-theme='c' class='list-info-delete'><a href='#'>"+i+": "+this.NOTE+"</a></li>";
							}
							else{
								dom = "<li data-icon='info' class='list-info'><a href='#'>"+i+": "+this.NOTE+"</a></li>";
							}
							$( "#" + info_dom_id ).append(dom).listview("refresh").fadeIn("fast");
						}
					});

					if( obj.SUBMISSION_STATUS == "success" ){
						$( ".highlight").removeClass("highlight");
						if ( action.toLowerCase().indexOf("download") >= 0){
							window.open(edge_path + obj.PATH);
						}else{
							if ( action.toLowerCase().indexOf("upload2gisaid") >= 0){
								obj.w = w;
								obj.spinner_id = 'newWindowSpinner';
								obj.type = "submit";
								obj.projname = projname;
								if (obj.PID) {
									checkpidInterval = setInterval(function(){check_process(obj)},3000);
								}

							}
							if (action.toLowerCase().indexOf("update") >= 0){
								$( "#edge_integrity_dialog_header" ).text("Message");
								$( "#edge_integrity_dialog_content" ).text("Your project(s) metadata was successfully updated.");
								setTimeout( function() { $( "#edge_integrity_dialog" ).popup('open'); }, 300 );
							}
							if (action.toLowerCase().indexOf("batch") >= 0 ){
								$( "#edge-gisaid-metadata-project-page" ).hide();
								$( "#edge-project-page" ).show();
							}else{
								setTimeout(function(){ updateReport($('#edge-output-projid').attr("data-pid"));},5000);
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
									var dom = "<li data-icon='delete' data-theme='c' class='list-info-delete'><a href='#'>"+v+"</a></li>";
									$( "#" + info_dom_id ).fadeIn("fast");
									$(dom).appendTo("#" + info_dom_id).on("click", function(){
										$('html, body').animate({
											scrollTop: $("#"+i).offset().top-100
										}, 200);
									});
								});
							});
						}else{
							var dom;
							if ( action === "upload2gisaid"){
								dom = "<li data-icon='delete' data-theme='c' class='list-info-delete'><a href='#'>FAILED to submit to GISAID & NCBI. Please check UPLOAD/submit.log in the project directory for detail.</a></li>";
							}
							if ( action === "batch-upload2gisaid"){
								dom = "<li data-icon='delete' data-theme='c' class='list-info-delete'><a href='#'>FAILED to submit to GISAID & NCBI.</a></li>";
							}
							$( "#" + info_dom_id ).append(dom).fadeIn("fast");
						}
						if (w){
							setTimeout(function(){ w.close(); },300);
						}
					}
				
					if( $(".list-info, .list-info-delete").size() ){
						var h = $( "#" + info_dom_id ).outerHeight();
						var h1 = parseInt(h)+100;
						$('html, body').animate({ scrollTop: "+="+h1+"px" }, 500);
					}

				},
				error: function(data){
					var dom = "<li data-icon='delete' data-theme='c' class='list-info-delete'><a href='#'>FAILED to run gisaid submission CGI. Please check server error log for detail.</a></li>";
					$( "#" +info_dom_id ).append(dom).fadeIn("fast");
					//listview("refresh");
					if (w){ w.close();}
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
				$( "#edge-gisaid-metadata-project-page" ).hide();		
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
	function check_process(data){
		var spinner_id = data.spinner_id;
		var w = data.w;
		$.ajax({
			url: "./cgi-bin/edge_action.cgi",
			type: "POST",
			dataType: "json",
			cache: false,
			data: { "action": 'checkpid', "pid": data.PID,'protocol': location.protocol, 'sid':localStorage.sid},
			success: function(obj){
				if (data.type === "submit"){
					var msgObj = $(w.document.body).find('#newWindowMsg');
					getLog('./' + data.PATH, msgObj);
				}
				if( obj.STATUS == "DONE" ){
					clearInterval(checkpidInterval);
					if (spinner_id){
						$('#' + spinner_id).removeClass("edge-sp edge-sp-circle");
					}else{
						$.mobile.loading( "hide");
					}
					if (data.type === "submit"){
						var spinnerObj = $(w.document.body).find('#newWindowSpinner');
						spinnerObj.removeClass("edge-sp edge-sp-circle");
						var msgObj = $(w.document.body).find('#newWindowMsg');
						var infoObj = $(w.document.body).find('#newWindowInfo');
						infoObj.text(function(){
							return $(this).text().replace("Please wait...","");
						});
						getLog('./' + data.PATH, msgObj, data.projname, true);
						//updateReport($('#edge-output-projid').attr("data-pid"));
 						//$('#get_download_link').after(data.LINK);
 						//$('#ddownload_link').addClass("ui-btn ui-mini ui-btn-inline ui-btn-active");
 					//	$('#ddownload_link').text('Download Project');
					//	$('#get_download_link').hide();
 					//	showWarning(data.INFO + data.LINK + '.');
					//	w.location = edge_path + data.PATH;
					//	setTimeout(function(){
					//		w.document.title = newTitle;
					//		var newHeader= "<div style='background:#02497e;'><h2 style='position:inherit; padding-left:20px;' class='edge-header'>"+newTitle+"</h2></div>";
					//		$(w.document.body).prepend(newHeader);
					//	},300);
					}else{
					//	w.opener.location = edge_path + data.PATH;
						setTimeout(function(){ w.close(); },300);
					}
					
				}else{
					//w.document.write(".");
					//console.log(obj.INFO);
				}
			},
			error: function(obj){
				showWarning("ACTION FAILED: Please try again or contact your system administrator.");
			}
		});
		
	}

	function getLog(texturl,dom,projname,end) {
                $.ajax({
                        url: texturl,
                        dataType: 'text',
                        cache: false,
                        success: function(text) {
                                dom.html(text);
				if (end){
					var msg;
					if ( /GISAID submit Completed/.test(text) && /NCBI submit Completed/.test(text)){
						msg = "Projects " + projname + " submmited to GISAID and NCBI sucessfully.";
					}else if ( /GISAID submit Completed/.test(text) ){
						msg = "Projects " + projname + " submmited to GISAID sucessfully. Please check log window for NCBI submission.";
					}else if ( /NCBI submit Completed/.test(text) ){
						if ( /sra_submit/.test(texturl) ){
							msg = "Projects " + projname + " submmited to NCBI SRA sucessfully.";
						}else{
							msg = "Projects " + projname + " submmited to NCBI sucessfully. Please check log window for GISAID submission.";
						}
					}else{
						if ( /sra_submit/.test(texturl) ){
							msg = "Projects " + projname + " submmited to NCBI SRA failed. Please check log window.";
						}else{
							msg = "Projects " + projname + " submission to GISAID and NCBI failed. Please check log window.";
						}
					}
					showWarning(msg);
				}
                        },
                        error: function(text) {
                                dom.text("No log retrieved from " + texturl);
                        }
                })
        };

	function showWarning( dialog_content ) {
		$( "#edge_integrity_dialog_content" ).html(dialog_content);
		$( "#edge_integrity_dialog" ).popup('open');
	}

	//NCBI SRA
	function check_bioproject_id (){
		if ( $("#metadata-sra-bioproject-sw1").is(':checked') ){
			$(".metadata-sra-bioproject-use-id").show();
			$(".metadata-sra-bioproject-input").hide();
		}
		if ( $("#metadata-sra-bioproject-sw2").is(':checked') ){
			$(".metadata-sra-bioproject-use-id").hide();
			$(".metadata-sra-bioproject-input").show();
		}
	}
	check_bioproject_id();
	$('#metadata-sra-bioproject-sw').on("change",check_bioproject_id);
	
	var model_optgroups = $('#metadata-sra-meta-libmodel').find('optgroup');
	$('#metadata-sra-meta-platform').on("change",function(){
		var platform = $(this).val();
		var previous_optgroups = $('#metadata-sra-meta-libmodel').find('optgroup').detach();
		var element = $.grep(model_optgroups,function(n,i){
					return n.label == platform;
			});
		$(element).appendTo($("#metadata-sra-meta-libmodel"));
        $('#metadata-sra-meta-libmodel').selectmenu("refresh");
	});
	$('#metadata-sra-meta-libmodel').on("change",function(){
		var platform = $("#metadata-sra-meta-libmodel option:selected").parent().attr('label');
		$('#metadata-sra-meta-platform').val(platform).selectmenu("refresh");
	});
	$( "#edge-sra-form-submit, #edge-sra-form-download, #edge-sra-form-update" ).on( "click", function() {
		var action = (this.id.toLowerCase().indexOf("download") >= 0)? "download" : (this.id.toLowerCase().indexOf("update") >= 0)? "update" : "upload2sra";
		var proj=$('#edge-output-projid').attr("data-pid");
		var projname=$('#edge-project-title').text().replace(' /','');
		var formDom = $("#edge-sra-upload-form");
		if ( action == 'upload2sra' ){
			var actionContent =  'By clicking the "Confirm" button, you hereby authorize EDGE-COVID19 to submit the samples and metadata to the NCBI SRA, and agree to remit the samples and related metadata to the public domain.</p>';
			$("#edge_confirm_dialog_content").html(actionContent);
			$('#edge_confirm_dialog').enhanceWithin().popup('open');
			$("#edge_confirm_dialog a:contains('Confirm')").unbind('click').on("click",function(){
				sra_actions(proj,formDom,action,projname);
            });
		}else{
			sra_actions(proj,formDom,action,projname);
		}
		
	});

	function sra_actions(proj, form, action, projname){
		var w;
		var info_dom_id = ( action.indexOf("batch") >=0 )? "edge-sra-batch-submit-info" : "edge-sra-submit-info";
		$.ajax({
				url: "./cgi-bin/edge_sra_upload.cgi",
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
					if (action.toLowerCase().indexOf("upload2sra") >= 0){
						w = window.open("","new","width=480,height=480");
						w.document.body.innerHTML = '';
						w.document.write( newWindowHeader + "<p id='newWindowInfo'>Submit project(s) " + projname + " to NCBI SRA. Please wait...</p>" + newWindowFooter);
					}
				},
				complete: function(data) {
					$("#" + info_dom_id).listview("refresh");
					setTimeout(function(){ $.mobile.loading( "hide" );},1000);
				},
				success: function(obj){
					// display general submission error
					$.each(obj, function(i,v){
						if( i!="PARAMS" && i!="SUBMISSION_STATUS" && i!="PATH" && i!="PID"){
							var dom;
							if( this.STATUS == "failure" ){
								dom = "<li data-icon='delete' data-theme='c' class='list-info-delete'><a href='#'>"+i+": "+this.NOTE+"</a></li>";
							}
							else{
								dom = "<li data-icon='info' class='list-info'><a href='#'>"+i+": "+this.NOTE+"</a></li>";
							}
							$( "#" + info_dom_id ).append(dom).listview("refresh").fadeIn("fast");
						}
					});

					if( obj.SUBMISSION_STATUS == "success" ){
						$( ".highlight").removeClass("highlight");
						if ( action.toLowerCase().indexOf("download") >= 0){
							window.open(edge_path + obj.PATH);
						}else{
							if ( action.toLowerCase().indexOf("upload2sra") >= 0){
								obj.w = w;
								obj.spinner_id = 'newWindowSpinner';
								obj.type = "submit";
								obj.projname = projname;
								if (obj.PID) {
									checkpidInterval = setInterval(function(){check_process(obj)},3000);
								}

							}
							if (action.toLowerCase().indexOf("update") >= 0){
								$( "#edge_integrity_dialog_header" ).text("Message");
								$( "#edge_integrity_dialog_content" ).text("Your project(s) SRA metadata was successfully updated.");
								setTimeout( function() { $( "#edge_integrity_dialog" ).popup('open'); }, 300 );
							}
							if (action.toLowerCase().indexOf("batch") >= 0 ){
								$( "#edge-sra-metadata-project-page" ).hide();
								$( "#edge-project-page" ).show();
							}else{
								setTimeout(function(){ updateReport($('#edge-output-projid').attr("data-pid"));},5000);
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
									var dom = "<li data-icon='delete' data-theme='c' class='list-info-delete'><a href='#'>"+v+"</a></li>";
									$( "#" + info_dom_id ).fadeIn("fast");
									$(dom).appendTo("#" + info_dom_id).on("click", function(){
										$('html, body').animate({
											scrollTop: $("#"+i).offset().top-100
										}, 200);
									});
								});
							});
						}else{
							var dom;
							if ( action === "upload2sra"){
								dom = "<li data-icon='delete' data-theme='c' class='list-info-delete'><a href='#'>FAILED to submit to NCBI SRA. Please check UPLOAD/sra_submit.log in the project directory for detail.</a></li>";
							}
							if ( action === "batch-upload2gisaid"){
								dom = "<li data-icon='delete' data-theme='c' class='list-info-delete'><a href='#'>FAILED to submit to NCBI SRA.</a></li>";
							}
							$( "#" + info_dom_id ).append(dom).fadeIn("fast");
						}
						if (w){
							setTimeout(function(){ w.close(); },300);
						}
					}
				
					if( $(".list-info, .list-info-delete").size() ){
						var h = $( "#" + info_dom_id ).outerHeight();
						var h1 = parseInt(h)+100;
						$('html, body').animate({ scrollTop: "+="+h1+"px" }, 500);
					}

				},
				error: function(data){
					var dom = "<li data-icon='delete' data-theme='c' class='list-info-delete'><a href='#'>FAILED to run NCBI SRA submission CGI. Please check server error log for detail.</a></li>";
					$( "#" +info_dom_id ).append(dom).fadeIn("fast");
					//listview("refresh");
					if (w){ w.close();}
				}
			});
	}
});
//# sourceURL=edge-output.js
