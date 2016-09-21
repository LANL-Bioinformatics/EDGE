$( document ).ready(function()
{
	var travels = 100;
	setTravels();
	setSymptoms();
	setSampleBlock();
	function setSampleBlock() {
		var sample_type = $('input[name= metadata-sample-type]:checked').val();
		if ( sample_type == 'human' ){
			$( "#metadata-host-human" ).fadeIn('fast');
			$( "#metadata-host-block" ).fadeIn('fast');
			$( "#metadata-host-h" ).fadeIn('fast');
			$( "#metadata-host-a" ).fadeOut('fast');

		}
		else if (sample_type == 'environmental' ){
			$( "#metadata-host-human" ).fadeOut('fast');
			$( "#metadata-host-block" ).fadeOut('fast');
			
		}
		else if (sample_type == 'animal' ){
			$( "#metadata-host-human" ).fadeOut('fast');
			$( "#metadata-host-block" ).fadeIn('fast');
			$( "#metadata-host-a" ).fadeIn('fast');
			$( "#metadata-host-h" ).fadeOut('fast');
		}
	}

	//checkbox
	if ($( '#metadata-host-gender-cb' ).is(':checked')) {
		$( "#human-gender" ).fadeIn('fast');
	} else {
		$( "#human-gender" ).fadeOut('fast');
	}
	$( '#metadata-host-gender-cb' ).on("change",function(){
		if($(this).is(':checked')){
			$( "#human-gender" ).fadeIn('fast');
		} else {
			$( "#human-gender" ).fadeOut('fast');
		}
	});
	if($( '#metadata-host-age-cb' ).is(':checked')){
		$( "#human-age" ).fadeIn('fast');
	} else {
		$( "#human-age" ).fadeOut('fast');
	}
	$( '#metadata-host-age-cb' ).on("change",function(){
		if($(this).is(':checked')){
			$( "#human-age" ).fadeIn('fast');
		} else {
			$( "#human-age" ).fadeOut('fast');
		}
	});

	//add a travel
   	$( '#add-travel' ).on("click",function(e){
        	e.preventDefault();
        	travels ++;
		var dom = '<div id="metadata-travel-'+travels+'">';
		dom += '<div class="ui-field-contain">';
		dom += '<label for="remove-travel">Travel #: </label>';
		dom += '<a href="#" class="ui-btn ui-btn-icon-left ui-icon-delete ui-mini" id="remove-travel">Remove this travel</a>';
		dom += '</div>';
								
		dom += '<div class="ui-field-contain">';
		dom += '<label>From</label>';
		dom += '<input data-role="date" type="text" data-mini="true" data-clear-btn="false" name="metadata-travel-date-f" id="metadata-travel-date-f-'+travels+'" maxlength="30" class="metadata-travel-date">';
		dom += '</div>';
		dom += '<div class="ui-field-contain">';
		dom += '<label>To</label>';
		dom += '<input data-role="date" type="text" data-mini="true" data-clear-btn="false" name="metadata-travel-date-t" id="metadata-travel-date-t-'+travels+'" maxlength="30" class="metadata-travel-date">';
		dom += '</div>';
		dom += '<div id="metadata-travel-geo-'+travels+'">';
		dom += '<div class="ui-field-contain">';
		dom += '<label>Location</label>';
		dom += '<input name="metadata-travel-location" id="geocomplete-travel-'+travels+'" type="text" placeholder="Type in an address to let system auto fill the location fields below"/>';
		dom += '</div>';
		dom += '<div class="ui-field-contain">';
		dom += '<label></label>';
		dom += '<input name="locality" id="metadata-travel-city-'+travels+'" data-mini="true" data-clear-btn="false" type="text" placeholder="City">';
		dom += '</div>';
		dom += '<div class="ui-field-contain">';
		dom += '<label></label>';
		dom += '<input name="administrative_area_level_1"  id="metadata-travel-state-'+travels+'" data-mini="true" data-clear-btn="false" type="text" placeholder="State">';
		dom += '</div>';
		dom += '<div class="ui-field-contain">';
		dom += '<label></label>';
		dom += '<input name="country" id="metadata-travel-country-'+travels+'" data-mini="true" data-clear-btn="false" type="text" placeholder="Country">';
		dom += '</div>';
		dom += '<div class="ui-field-contain">';
		dom += '<label></label>';
		dom += '<input name="lat" id="metadata-travel-lat-'+travels+'" data-mini="true" data-clear-btn="false" type="text" placeholder="Latitude">';
		dom += '</div>';
		dom += '<div class="ui-field-contain">';
		dom += '<label></label>';
		dom += '<input name="lng" id="metadata-travel-lng-'+travels+'" data-mini="true" data-clear-btn="false" type="text" placeholder="Longitude">';
		dom += '</div>';
		dom += '</div>';
		dom += '<br><br>';
		dom += '</div>';
		$("#metadata-travels" ).append(dom).trigger("create");
	});

	//remove a travel
	$( '#metadata-travels' ).on("click","#remove-travel",function(){
		$(this).parent('div').parent('div').remove();
	});
	//travel geo location, limit to 10
	$('#metadata-travels').on('click','#geocomplete-travel-1',function() {
		$(this).geocomplete({
		  details: "#metadata-travels #metadata-travel-geo-1"
		});
	});
	$('#metadata-travels').on('click','#geocomplete-travel-2',function() {
		$(this).geocomplete({
		  details: "#metadata-travels #metadata-travel-geo-2"
		});
	});
	$('#metadata-travels').on('click','#geocomplete-travel-3',function() {
		$(this).geocomplete({
		  details: "#metadata-travels #metadata-travel-geo-3"
		});
	});
	$('#metadata-travels').on('click','#geocomplete-travel-4',function() {
		$(this).geocomplete({
		  details: "#metadata-travels #metadata-travel-geo-4"
		});
	});
	$('#metadata-travels').on('click','#geocomplete-travel-5',function() {
		$(this).geocomplete({
		  details: "#metadata-travels #metadata-travel-geo-5"
		});
	});
	$('#metadata-travels').on('click','#geocomplete-travel-6',function() {
		$(this).geocomplete({
		  details: "#metadata-travels #metadata-travel-geo-6"
		});
	});
	$('#metadata-travels').on('click','#geocomplete-travel-7',function() {
		$(this).geocomplete({
		  details: "#metadata-travels #metadata-travel-geo-7"
		});
	});
	$('#metadata-travels').on('click','#geocomplete-travel-8',function() {
		$(this).geocomplete({
		  details: "#metadata-travels #metadata-travel-geo-8"
		});
	});
	$('#metadata-travels').on('click','#geocomplete-travel-9',function() {
		$(this).geocomplete({
		  details: "#metadata-travels #metadata-travel-geo-9"
		});
	});
	$('#metadata-travels').on('click','#geocomplete-travel-10',function() {
		$(this).geocomplete({
		  details: "#metadata-travels #metadata-travel-geo-10"
		});
	});

	$('#metadata-travels').on('click','#geocomplete-travel-101',function() {
		$(this).geocomplete({
		  details: "#metadata-travels #metadata-travel-geo-101"
		});
	});
	$('#metadata-travels').on('click','#geocomplete-travel-102',function() {
		$(this).geocomplete({
		  details: "#metadata-travels #metadata-travel-geo-102"
		});
	});
	$('#metadata-travels').on('click','#geocomplete-travel-103',function() {
		$(this).geocomplete({
		  details: "#metadata-travels #metadata-travel-geo-103"
		});
	});
	$('#metadata-travels').on('click','#geocomplete-travel-4',function() {
		$(this).geocomplete({
		  details: "#metadata-travels #metadata-travel-geo-104"
		});
	});
	$('#metadata-travels').on('click','#geocomplete-travel-105',function() {
		$(this).geocomplete({
		  details: "#metadata-travels #metadata-travel-geo-105"
		});
	});
	$('#metadata-travels').on('click','#geocomplete-travel-106',function() {
		$(this).geocomplete({
		  details: "#metadata-travels #metadata-travel-geo-106"
		});
	});
	$('#metadata-travels').on('click','#geocomplete-travel-107',function() {
		$(this).geocomplete({
		  details: "#metadata-travels #metadata-travel-geo-107"
		});
	});
	$('#metadata-travels').on('click','#geocomplete-travel-108',function() {
		$(this).geocomplete({
		  details: "#metadata-travels #metadata-travel-geo-108"
		});
	});
	$('#metadata-travels').on('click','#geocomplete-travel-109',function() {
		$(this).geocomplete({
		  details: "#metadata-travels #metadata-travel-geo-109"
		});
	});
	$('#metadata-travels').on('click','#geocomplete-travel-110',function() {
		$(this).geocomplete({
		  details: "#metadata-travels #metadata-travel-geo-110"
		});
	});
	$('#metadata-travels').on('focus','.metadata-travel-date',function() {
		$(this).datepicker({
		      changeMonth: true,
		      changeYear: true,
		      dateFormat: 'yy-mm-dd'
	    	});
	});

	//sample type change
	$('input[name= metadata-sample-type]' ).on("change",function(){
		$( "#metadata-isolation-source" ).val("");  
		if ( $(this).val() == 'human' ){
			$( "#metadata-host-human" ).fadeIn('fast');
			$( "#metadata-host-block" ).fadeIn('fast');
			$( "#metadata-host-h" ).fadeIn('fast');
			$( "#metadata-host-a" ).fadeOut('fast');

		}
		else if ($(this).val() == 'environmental' ){
			$( "#metadata-host-human" ).fadeOut('fast');
			$( "#metadata-host-block" ).fadeOut('fast');
			
		}
		else if ($(this).val() == 'animal' ){
			$( "#metadata-host-human" ).fadeOut('fast');
			$( "#metadata-host-block" ).fadeIn('fast');
			$( "#metadata-host-a" ).fadeIn('fast');
			$( "#metadata-host-h" ).fadeOut('fast');
		}
	});	

	//geo location
	$('#metadata-sample-geo').on('click','#metadata-sample-geocomplete',function() {
		$(this).geocomplete({
		  details: "#metadata-sample-geo"
		});
	});

	//datepicker
	$('.metadata-input-date').datepicker({
		changeMonth: true,
		changeYear: true,
		dateFormat: 'yy-mm-dd'
	});

	//study-title auto complete
	$(".dblist-study ul").on('click', 'li', function() {
		$( "#metadata-study-title" ).val( $( this ).text() );  
		$("#metadata-study-title-list li" ).addClass('ui-screen-hidden');
		$('#study-title-div').css('display', 'none');
	});

	$("#metadata-study-title").focus(function(){
		setStudyList();
		$('#study-title-div').css('display', 'inline');
	});

	$("#metadata-study-title").focusout(function(){
		if ($("#metadata-study-title-list li:hover").length == 0) {
			$("#metadata-study-title-list li" ).addClass('ui-screen-hidden');
			$('#study-title-div').css('display', 'none');
		}
	});

	//study-type auto complete
	$(".dblist-study-type ul").on('click', 'li', function() {
		$( "#metadata-study-type" ).val( $( this ).text() );  
		$("#metadata-study-type-list li" ).addClass('ui-screen-hidden');
		$('#study-type-div').css('display', 'none');
	});

	$("#metadata-study-type").focus(function(){
		setStudyTypeList();
		$('#study-type-div').css('display', 'inline');
	});

	$("#metadata-study-type").focusout(function(){
		if ($("#metadata-study-type-list li:hover").length == 0) {
			$("#metadata-study-type-list li" ).addClass('ui-screen-hidden');
			$('#study-type-div').css('display', 'none');
		}
	});

	//pg host auto complete
	$(".dblist-pg-host ul").on('click', 'li', function() {
		$( "#metadata-host" ).val( $( this ).text() );  
		$("#metadata-host-list li" ).addClass('ui-screen-hidden');
		$('#host-list-div').css('display', 'none');
	});

	$("#metadata-host").focus(function(){
		setPgHostList();
		$('#host-list-div').css('display', 'inline');
	});

	$("#metadata-host").focusout(function(){
		if ($("#metadata-host-list li:hover").length == 0) {
			$("#metadata-host-list li" ).addClass('ui-screen-hidden');
			$('#host-list-div').css('display', 'none');
		}
	});

	//isolation source auto complete
	$(".dblist-sample-source ul").on('click', 'li', function() {
		$( "#metadata-isolation-source" ).val( $( this ).text() );  
		$("#metadata-isolation-source-list li" ).addClass('ui-screen-hidden');
		$('#isolation-source-list-div').css('display', 'none');	
	});

	$("#metadata-isolation-source").focus(function(){
		setIsolationSourceList($('input[name= metadata-sample-type]:checked').val());
		$('#isolation-source-list-div').css('display', 'inline');		
	});

	$("#metadata-isolation-source").focusout(function(){
		if ($("#metadata-isolation-source-list li:hover").length == 0) {
			$("#metadata-isolation-source-list li" ).addClass('ui-screen-hidden');
			$('#isolation-source-list-div').css('display', 'none');	
		}
	});

	//sequencing center auto complete
	$(".dblist-seq-center ul").on('click', 'li', function() {
		$( "#metadata-seq-center" ).val( $( this ).text() );  
		$("#edge-seq-center-list li" ).addClass('ui-screen-hidden');
		$('#seq-center-list-div').css('display', 'none');
	});

	$( "#metadata-seq-center" ).focus(function(){
		setSeqCenterList();
		$('#seq-center-list-div').css('display', 'inline');
	});

	$( "#metadata-seq-center" ).focusout(function(){
		if ($("#edge-seq-center-list li:hover").length == 0) {
			$("#edge-seq-center-list li" ).addClass('ui-screen-hidden');
			$('#seq-center-list-div').css('display', 'none');
		}
	});

	//sequencer auto complete
	$(".dblist-sequencer ul").on('click', 'li', function() {
		$( "#metadata-sequencer" ).val( $( this ).text() );  
		$("#edge-sequencer-list li" ).addClass('ui-screen-hidden');
		$('#sequencer-list-div').css('display', 'none');
	});

	$( "#metadata-sequencer" ).focus(function(){
		setSequencerList();
		$('#sequencer-list-div').css('display', 'inline');
	});

	$( "#metadata-sequencer" ).focusout(function(){
		if ($("#edge-sequencer-list li:hover").length == 0) {
			$("#edge-sequencer-list li" ).addClass('ui-screen-hidden');
			$('#sequencer-list-div').css('display', 'none');
		}
	});

	function setStudyList(){
    		$.ajax({
    			type: "POST",
			url: "./cgi-bin/edge_db.cgi", 
        		cache: false,
        		dataType: "html",
			data: {"action":"study-list" },
        		// script call was successful 
        		success: function(data){
				$("#metadata-study-title-list" ).empty();
				var its = data.split("\n");
				for (var i=0;i<its.length;i++) {
					$("#metadata-study-title-list" ).append($("<li>"+its[i]+"</it>"));
				}
				$("#metadata-study-title-list li" ).addClass('ui-screen-hidden');
       			} // success
    		});
	}

	function setStudyTypeList(){
    		$.ajax({
    			type: "POST",
			url: "./cgi-bin/edge_db.cgi", 
        		cache: false,
        		dataType: "html",
			data: {"action":"study-type-list" },
        		// script call was successful 
        		success: function(data){
				$("#metadata-study-type-list" ).empty();
				var its = data.split("\n");
				for (var i=0;i<its.length;i++) {
					$("#metadata-study-type-list" ).append($("<li>"+its[i]+"</it>"));
				}
				$("#metadata-study-type-list li" ).addClass('ui-screen-hidden');
       			} // success
    		});
	}

	function setPgHostList(){
    		$.ajax({
    			type: "POST",
			url: "./cgi-bin/edge_db.cgi", 
        		cache: false,
        		dataType: "html",
			data: {"action":"pg-host-list" },
        		// script call was successful 
        		success: function(data){
				$("#metadata-host-list" ).empty();
				var its = data.split("\n");
				for (var i=0;i<its.length;i++) {
					$("#metadata-host-list" ).append($("<li>"+its[i]+"</it>"));
				}
				$("#metadata-host-list li" ).addClass('ui-screen-hidden');
       			} // success
    		});
	}

	function setIsolationSourceList(type){
    		$.ajax({
    			type: "POST",
			url: "./cgi-bin/edge_db.cgi", 
        		cache: false,
        		dataType: "html",
			data: {"action":"isolation-source-list", "sample-type":type },
        		// script call was successful 
        		success: function(data){
				$("#metadata-isolation-source-list" ).empty();
				var its = data.split("\n");
				for (var i=0;i<its.length;i++) {
					$("#metadata-isolation-source-list" ).append($("<li>"+its[i]+"</it>"));
				}
				$("#metadata-isolation-source-list li" ).addClass('ui-screen-hidden');
       			} // success
    		});
	}

	function setSeqCenterList(){
    		$.ajax({
    			type: "POST",
			url: "./cgi-bin/edge_db.cgi", 
        		cache: false,
        		dataType: "html",
			data: {"action":"seq-center-list" },
        		// script call was successful 
        		success: function(data){
				$("#edge-seq-center-list" ).empty();
				var its = data.split("\n");
				for (var i=0;i<its.length;i++) {
					$("#edge-seq-center-list" ).append($("<li>"+its[i]+"</it>"));
				}
				$("#edge-seq-center-list li" ).addClass('ui-screen-hidden');
       			} // success
    		});
	}

	function setSequencerList(){
    		$.ajax({
    			type: "POST",
			url: "./cgi-bin/edge_db.cgi", 
        		cache: false,
        		dataType: "html",
			data: {"action":"sequencer-list" },
        		// script call was successful 
        		success: function(data){
				$("#edge-sequencer-list" ).empty();
				var its = data.split("\n");
				for (var i=0;i<its.length;i++) {
					$("#edge-sequencer-list" ).append($("<li>"+its[i]+"</it>"));
				}
				$("#edge-sequencer-list li" ).addClass('ui-screen-hidden');
       			} // success
    		});
	}

	function setTravels(){
		var proj=$('#edge-output-projid').attr("data-pid")
    		$.ajax({
    			type: "POST",
			url: "./cgi-bin/edge_db.cgi", 
        		cache: false,
        		dataType: "html",
			data: {"action":"travel-edit-list", "proj": proj},
        		// script call was successful 
        		success: function(data){
				$("#metadata-travels" ).empty();
				$("#metadata-travels" ).append(data).trigger("create");
       			}, // success
			error: function(data){
				console.log(data);
			}
    		});
	}

	function setSymptoms(){
		var proj=$('#edge-output-projid').attr("data-pid")
    		$.ajax({
    			type: "POST",
			url: "./cgi-bin/edge_db.cgi", 
        		cache: false,
        		dataType: "html",
			data: {"action":"symptom-edit-list", "proj": proj},
        		// script call was successful 
        		success: function(data){
				$("#metadata-symptoms" ).empty();
				$("#metadata-symptoms" ).append(data).trigger("create");
       			}, // success
			error: function(data){
				console.log(data);
			}
    		});
	}


	//form reset
	$( "#edge-sample-metadata-form-reset" ).on( "click", function() {
		$( "#edge-sample-metadata-form-edit" )[0].reset();

		//sample type
		setSampleBlock();
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
});
//# sourceURL=edge-output.js
