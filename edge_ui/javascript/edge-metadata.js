$( document ).ready(function()
{
	//gender and age
	if($('#pg-cb-gender-edit').is(':checked')){
		$( "#human-gender-edit" ).fadeIn('fast');
	} else {
		$( "#human-gender-edit" ).fadeOut('fast');
	}
	if($('#pg-cb-age-edit').is(':checked')){
		$( "#human-age-edit" ).fadeIn('fast');
	} else {
		$( "#human-age-edit" ).fadeOut('fast');
	}

	$( '#pg-cb-gender-edit' ).on("change",function(){
		if($(this).is(':checked')){
			$( "#human-gender-edit" ).fadeIn('fast');
		} else {
			$( "#human-gender-edit" ).fadeOut('fast');
		}
	});
	$( '#pg-cb-age-edit' ).on("change",function(){
		if($(this).is(':checked')){
			$( "#human-age-edit" ).fadeIn('fast');
		} else {
			$( "#human-age-edit" ).fadeOut('fast');
		}
	});

	//sample type
	if($( '#edge-sample-type-edit').val() == "human") {
			$( "#edge-pg-human-edit" ).fadeIn('fast');
			$( "#edge-pg-host-associated-edit" ).fadeIn('fast');
			$( "#edge-sample-source-options-host-edit" ).fadeIn('fast');
			$( "#edge-sample-source-options-nonhost-edit" ).fadeOut('fast');
	}
	if($( '#edge-sample-type-edit' ).val() == "animal") {
			$( "#edge-pg-human-edit" ).fadeOut('fast');
			$( "#edge-pg-host-associated-edit" ).fadeIn('fast');
			$( "#edge-sample-source-options-host-edit" ).fadeIn('fast');
			$( "#edge-sample-source-options-nonhost-edit" ).fadeOut('fast');
	}
	if($( '#edge-sample-type-edit' ).val() == "environmental") {
			$( "#edge-pg-human-edit" ).fadeOut('fast');
			$( "#edge-pg-host-associated-edit" ).fadeOut('fast');
			$( "#edge-sample-source-options-host-edit" ).fadeOut('fast');
			$( "#edge-sample-source-options-nonhost-edit" ).fadeIn('fast');
	}

	$( '#edge-sample-type-edit' ).on("change",function(){
		if ( $(this).val() == 'human' ){
			$( "#edge-pg-human-edit" ).fadeIn('fast');
			$( "#edge-pg-host-associated-edit" ).fadeIn('fast');
			$('#edge-pg-host-edit').val("human");
			$( "#edge-sample-source-options-host-edit" ).fadeIn('fast');
			$( "#edge-sample-source-options-nonhost-edit" ).fadeOut('fast');

		}
		else if ($(this).val() == 'environmental' ){
			$( "#edge-pg-human-edit" ).fadeOut('fast');
			$( "#edge-pg-host-associated-edit" ).fadeOut('fast');
			$( "#edge-sample-source-options-host-edit" ).fadeOut('fast');
			$( "#edge-sample-source-options-nonhost-edit" ).fadeIn('fast');
			
		}
		else if ($(this).val() == 'animal' ){
			$( "#edge-pg-human-edit" ).fadeOut('fast');
			$( "#edge-pg-host-associated-edit" ).fadeIn('fast');
			$('#edge-pg-host-edit').val("");
			$( "#edge-sample-source-options-host-edit" ).fadeIn('fast');
			$( "#edge-sample-source-options-nonhost-edit" ).fadeOut('fast');
		}
	});

	//seq platform
	if($( '#edge-pg-seq-platform-edit' ).val() =='Illumina') {
			$( "#edge-pg-sequencer-options-ill-edit" ).fadeIn('fast');
			$( "#edge-pg-sequencer-options-ion-edit" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-nan-edit" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-pac-edit" ).fadeOut('fast');

	}
	if($( '#edge-pg-seq-platform-edit' ).val() =='IonTorrent') {
			$( "#edge-pg-sequencer-options-ill-edit" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-ion-edit" ).fadeIn('fast');
			$( "#edge-pg-sequencer-options-nan-edit" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-pac-edit" ).fadeOut('fast');

	}
	if($( '#edge-pg-seq-platform-edit' ).val() =='Nanopore') {
			$( "#edge-pg-sequencer-options-ill-edit" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-ion-edit" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-nan-edit" ).fadeIn('fast');
			$( "#edge-pg-sequencer-options-pac-edit" ).fadeOut('fast');

	}
	if($( '#edge-pg-seq-platform-edit' ).val() =='PacBio') {
			$( "#edge-pg-sequencer-options-ill-edit" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-ion-edit" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-nan-edit" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-pac-edit" ).fadeIn('fast');

	}

	$( '#edge-pg-seq-platform-edit' ).on("change",function(){
		if ( $(this).val() == 'Illumina' ){
			$( "#edge-pg-sequencer-options-ill-edit" ).fadeIn('fast');
			$( "#edge-pg-sequencer-options-ion-edit" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-nan-edit" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-pac-edit" ).fadeOut('fast');
		}
		else if ($(this).val() == 'IonTorrent' ){
			$( "#edge-pg-sequencer-options-ill-edit" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-ion-edit" ).fadeIn('fast');
			$( "#edge-pg-sequencer-options-nan-edit" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-pac-edit" ).fadeOut('fast');
		}
		else if ($(this).val() == 'Nanopore' ){
			$( "#edge-pg-sequencer-options-ill-edit" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-ion-edit" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-nan-edit" ).fadeIn('fast');
			$( "#edge-pg-sequencer-options-pac-edit" ).fadeOut('fast');
		}
		else if ($(this).val() == 'PacBio' ){
			$( "#edge-pg-sequencer-options-ill-edit" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-ion-edit" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-nan-edit" ).fadeOut('fast');
			$( "#edge-pg-sequencer-options-pac-edit" ).fadeIn('fast');
		}
	});

	//form reset
	$( "#edge-sample-metadata-form-reset" ).on( "click", function() {
		$( "#edge-sample-metadata-form-edit" )[0].reset();

		//sample type
		if($( '#edge-sample-type-edit').val() == "human") {
				$( "#edge-pg-human-edit" ).fadeIn('fast');
				$( "#edge-pg-host-associated-edit" ).fadeIn('fast');
				$( "#edge-sample-source-options-host-edit" ).fadeIn('fast');
				$( "#edge-sample-source-options-nonhost-edit" ).fadeOut('fast');
		}
		if($( '#edge-sample-type-edit' ).val() == "animal") {
				$( "#edge-pg-human-edit" ).fadeOut('fast');
				$( "#edge-pg-host-associated-edit" ).fadeIn('fast');
				$( "#edge-sample-source-options-host-edit" ).fadeIn('fast');
				$( "#edge-sample-source-options-nonhost-edit" ).fadeOut('fast');
		}
		if($( '#edge-sample-type-edit' ).val() == "environmental") {
				$( "#edge-pg-human-edit" ).fadeOut('fast');
				$( "#edge-pg-host-associated-edit" ).fadeOut('fast');
				$( "#edge-sample-source-options-host-edit" ).fadeOut('fast');
				$( "#edge-sample-source-options-nonhost-edit" ).fadeIn('fast');
		}
	});

	//form cancel
	$( "#edge-sample-metadata-form-cancel" ).on( "click", function() {
		updateReport($('#edge-output-projid').attr("data-pid"));
	});


	//form submit
	$( "#edge-sample-metadata-form-submit" ).on( "click", function() {

		var proj=$('#edge-output-projid').attr("data-pid")
		$.ajax({
			url: "./cgi-bin/edge_sample_metadata.cgi",
			type: "POST",
			dataType: "json",
			cache: false,
			//data: $( "#edge-run-pipeline-form" ).serialize(),
			data: ( $("#edge-sample-metadata-form-edit").serialize() +'&'+ $.param({ 'proj': proj, 'sid':localStorage.sid, 'action':'save' })),
			beforeSend: function(){
				$.mobile.loading( "show", {
					text: "submitting...",
					textVisible: 1,
					html: ""
				});
			},
			complete: function(data) {
				$.mobile.loading( "hide" );
				updateReport($('#edge-output-projid').attr("data-pid"));
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
					updateReport($('#edge-output-projid').attr("data-pid"));
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

	//geo location
        $("#geocomplete").geocomplete({
          map: ".map_canvas",
          details: "form",
          types: ["geocode", "establishment"],
        });
});
//# sourceURL=edge-output.js
