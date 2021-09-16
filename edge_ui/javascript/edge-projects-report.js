$( document ).ready(function()
{
	//checkbox
	//metadata
	if ($( '#checkbox-report-1-4' ).is(':checked')) {
		$( "#report-metadata-list" ).fadeIn('fast');
	} else {
		$( "#report-metadata-list" ).fadeOut('fast');
	}

	$( '#checkbox-report-1-4' ).on("change",function(){
		if($(this).is(':checked')){
			$( "#report-metadata-list" ).fadeIn('fast');
		} else {
			$( "#report-metadata-list" ).fadeOut('fast');
		}
	});
	//pre-processing 
	if ($( '#checkbox-report-2-1' ).is(':checked')) {
		$( "#report-preprocess-stats-list" ).fadeIn('fast');
	} else {
		$( "#report-preprocess-stats-list" ).fadeOut('fast');
	}

	$( '#checkbox-report-2-1' ).on("change",function(){
		if($(this).is(':checked')){
			$( "#report-preprocess-stats-list" ).fadeIn('fast');
		} else {
			$( "#report-preprocess-stats-list" ).fadeOut('fast');
		}
	});
	
	if ($( '#checkbox-report-2-2' ).is(':checked')) {
		$( "#report-preprocess-figure-list" ).fadeIn('fast');
	} else {
		$( "#report-preprocess-figure-list" ).fadeOut('fast');
	}

	$( '#checkbox-report-2-2' ).on("change",function(){
		if($(this).is(':checked')){
			$( "#report-preprocess-figure-list" ).fadeIn('fast');
		} else {
			$( "#report-preprocess-figure-list" ).fadeOut('fast');
		}
	});
	//assembly and annotation
	if ($( '#checkbox-report-3-1' ).is(':checked')) {
		$( "#report-aa-stats-list" ).fadeIn('fast');
	} else {
		$( "#report-aa-stats-list" ).fadeOut('fast');
	}

	$( '#checkbox-report-3-1' ).on("change",function(){
		if($(this).is(':checked')){
			$( "#report-aa-stats-list" ).fadeIn('fast');
		} else {
			$( "#report-aa-stats-list" ).fadeOut('fast');
		}
	});
	
	if ($( '#checkbox-report-3-2' ).is(':checked')) {
		$( "#report-aa-figure-list" ).fadeIn('fast');
	} else {
		$( "#report-aa-figure-list" ).fadeOut('fast');
	}

	$( '#checkbox-report-3-2' ).on("change",function(){
		if($(this).is(':checked')){
			$( "#report-aa-figure-list" ).fadeIn('fast');
		} else {
			$( "#report-aa-figure-list" ).fadeOut('fast');
		}
	});
	//Reference-Based Analysis
	if ($( '#checkbox-report-4-1' ).is(':checked')) {
		$( "#report-ref-stats-list" ).fadeIn('fast');
	} else {
		$( "#report-ref-stats-list" ).fadeOut('fast');
	}

	$( '#checkbox-report-4-1' ).on("change",function(){
		if($(this).is(':checked')){
			$( "#report-ref-stats-list" ).fadeIn('fast');
		} else {
			$( "#report-ref-stats-list" ).fadeOut('fast');
		}
	});
	
	if ($( '#checkbox-report-4-2' ).is(':checked')) {
		$( "#report-ref-figure-list" ).fadeIn('fast');
	} else {
		$( "#report-ref-figure-list" ).fadeOut('fast');
	}

	$( '#checkbox-report-4-2' ).on("change",function(){
		if($(this).is(':checked')){
			$( "#report-ref-figure-list" ).fadeIn('fast');
		} else {
			$( "#report-ref-figure-list" ).fadeOut('fast');
		}
	});
	//Taxonomy Classification read-based
	if ($( '#checkbox-report-5-1' ).is(':checked')) {
		$( "#report-tax-stats-list" ).fadeIn('fast');
	} else {
		$( "#report-tax-stats-list" ).fadeOut('fast');
	}

	$( '#checkbox-report-5-1' ).on("change",function(){
		if($(this).is(':checked')){
			$( "#report-tax-stats-list" ).fadeIn('fast');
		} else {
			$( "#report-tax-stats-list" ).fadeOut('fast');
		}
	});
	
	if ($( '#checkbox-report-5-2' ).is(':checked')) {
		$( "#report-tax-figure-list" ).fadeIn('fast');
	} else {
		$( "#report-tax-figure-list" ).fadeOut('fast');
	}

	$( '#checkbox-report-5-2' ).on("change",function(){
		if($(this).is(':checked')){
			$( "#report-tax-figure-list" ).fadeIn('fast');
		} else {
			$( "#report-tax-figure-list" ).fadeOut('fast');
		}
	});
	
	if ($( '#checkbox-report-5-3' ).is(':checked')) {
		$( "#report-tax-tool-list" ).fadeIn('fast');
	} else {
		$( "#report-tax-tool-list" ).fadeOut('fast');
	}

	$( '#checkbox-report-5-3' ).on("change",function(){
		if($(this).is(':checked')){
			$( "#report-tax-tool-list" ).fadeIn('fast');
		} else {
			$( "#report-tax-tool-list" ).fadeOut('fast');
		}
	});

	//form reset - not implemented yet
	$( "#edge-projects-report-form-reset" ).on( "click", function() {
		$( "#edge-report-form" )[0].reset();
		//need reset all controlgroup checkboxes and hide them
		//need unselect projects
	});

	//form cancel
	$( "#edge-projects-report-form-cancel" ).on( "click", function() {
		location.reload();
	});

	//form submit
	$( "#edge-projects-report-form-submit" ).on( "click", function() {
		var report=$('#edge-report-name').val();
		report = report.replace(/\s+/g, '');
		if(report==='') {
			showWarning("Report name is required.");
			return;
		}

		if ( $('[name="edge-reportform-ckb"]:checked').length === 0 ){
			showWarning("There are no projects selected.");
			return;
		}

		$.ajax({
			url: "./cgi-bin/edge_projects_report.cgi",
			type: "POST",
			dataType: "html",
			cache: false,
			data: ( $("#edge-report-form").serialize() +'&'+ $.param({ 'report': report, 'sid':localStorage.sid, 'action':'create' })),
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
				//console.log(data);
				var allMainPage = $(".edge-main-page");
				allMainPage.hide();
				$( "#edge-projects-report-page" ).html(data);
				$( "#edge-projects-report-page" ).show();

				$(".tooltip").tooltipster({
					theme:'tooltipster-light',
					maxWidth: '480',
					interactive: true,
					multiple:true
				});
				$( "#edge-projects-report-page" ).enhanceWithin();

				$.getScript( "./javascript/edge-projects-report.js" )
					.done(function( script, textStatus ) {
					//	console.log( "edge-output.js loaded: " + textStatus );
					})
					.fail(function( jqxhr, settings, exception ) {
						console.log( jqxhr, settings, exception );
					});
				//check report plot
				var reportURL = "." + $('#report-dir').text() + "/lanl_project_list_ec19.html";
				var reportLogURL = "." + $('#report-dir').text() + "/ec19_report.log";
				if (/EDGE_report/.test(location.pathname.toString()) || ! location.hostname){
					reportURL = "./lanl_project_list_ec19.html";
				}
				$('#report-plots-iframe').hide();
				$('#report-plots-link').hide();
				$('#wait-report-plots').show();
				var updatePlotInterval = setInterval(function(){
				$.get(reportURL)
					.done(function() { 
						clearInterval(updatePlotInterval);
						$('#report-plots-iframe').attr("src",reportURL).show();
						$('#report-plots-link').attr("href",reportURL).show();
						$('#wait-report-plots').hide();
					}).fail(function() { 
						$('#report-plots-iframe').hide();
						$('#report-plots-link').hide();
						$('#wait-report-plots').show();
						$.get(reportLogURL).done(function(data){
								if (/Error/.test(data)){
										clearInterval(updatePlotInterval);
										$( "#edge_integrity_dialog_content" ).text("Failed to generate the ec19 variant vis report. Please check the ec19_report.log");
										$( "#edge_integrity_dialog" ).popup('open');
								}
						});
					})
				}, 5000);
			},
			error: function(data){
				//console.log(data);
				$.mobile.loading( "hide" );
				$( "#edge_integrity_dialog_content" ).text("Failed to generate the report. Please REFRESH the page and try again.");
				$( "#edge_integrity_dialog" ).popup('open');
			}
		});
	});

	function showWarning( dialog_content ) {
		$( "#edge_integrity_dialog_content" ).html(dialog_content);
		$( "#edge_integrity_dialog" ).popup('open');
	}


});
