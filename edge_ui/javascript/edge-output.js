$( document ).ready(function()
{
	var checkpidInterval;
	var projdir = $('#edge-proj-outdir').attr('dir-src');
	var projNotes = projdir + "/projnotes.txt";
	var focusProjName = $('#edge-output-projname').attr("data-pid");
	var projRealname = $('#edge-output-projname').text();
	var newWindowHeader = "<html><head><title>EDGE bioinformatics</title><link rel='stylesheet' href='css/edge-output.css'/></head><div style='background:#50a253;'><h2 style='position:inherit; padding-left:20px;'>EDGE bioinformatics</h2></div>";
        var newWindowFooter = "<div class='edge-sp edge-sp-circle'></div></body></html>";
	var newTitle = "EDGE bioinformatics";


	$( "#edge-out-expand-all" ).on("click", function(){
		$( "#edge-content-report > div[data-role='collapsible']" ).collapsible( "option", "collapsed", false );
	});

	$( "#edge-out-expand-none" ).on("click", function(){
		$( "#edge-content-report > div[data-role='collapsible']" ).collapsible( "option", "collapsed", true );
	});
	$( ".edge-sgout-expand-all" ).on("click", function(){
		$(this).parent().next().find('.li-report-content').toggle(true);
	});
	$( ".edge-sgout-expand-none" ).on("click", function(){
		$(this).parent().next().find('.li-report-content').toggle(false);
	});
	$( "#edge-content-report p.extra-text-link > a:not(a[data-rel='popup'])" ).attr( "target", "new_window");
	$( "#edge-content-report .iframe_label > a:not(a[data-rel='popup'])" ).attr( "target", "new_window");
	$( "#edge-content-report .jb_link" ).attr( "target", "new_window");
	$( "#edge-content-report table.output-table td > a:not(a[data-rel='popup'],.anchorlink)" ).attr( "target", "new_window");
	$( ".edge-compare-output .iframe_label > a:not(a[data-rel='popup'])" ).attr( "target", "_blank");

	$( "#edge-content-report span.li-report-content-title" ).on("click",function(){
		$(this).next().slideToggle();
	});
	$( ".anchorlink").on("click",function(e){
		e.preventDefault();
		$('html, body').animate({
			scrollTop: $($(this).attr('href')).offset().top - 80 }, 300, 'linear'
		);
	});
	$( "#edge-content-report a.open-jbrowse-iframe" ).on("click",function(){
		var iobg = $(this).next();
		if( iobg.prop("tagName") == "IFRAME" && iobg.is(":visible") ){
			$(this).next().hide();
		}
		else if ( iobg.prop("tagName") == "IFRAME" ){
			$(this).next().show();
		}
		else{
			var url = $(this).prev().attr("href");
			var jb = $( '<iframe style="height:540px;"></iframe>' ).attr("src",url);
			$(this).after(jb);
		}
	});
	$( "#edge-content-report a.iframe_data_src" ).on("click",function(){
		var src = $(this).attr("data-src");
		$(this).siblings(':first').attr("href",src);
		$(this).parent().next('iframe').attr("src",src);
	});
	
	// directory file tree
	var dir;
	var loc = window.location.pathname.replace("//","/");
	var edge_path = loc.substring(0,loc.lastIndexOf('/'));
	$('#edge-outputfile-dialog').popup({ positionTo: "window"});
	$('.edge-outputfile-tree').on('click',function(){
		dir=$(this).attr('dir-src');
		$('#edge-outputfile-tree').fileTree({
			root: '/' + dir + '/',
			script: './cgi-bin/jqueryFileTree_rel.cgi',
			}, function(file) {
				window.open(edge_path + file);
		});
	});
	//data Table
	 $('.SGtable').removeClass("ui-table-reflow ui-table");
	 $('.SGtable').find(".ui-table-cell-label").remove();
	 $('.SGtable').each(function(){
                var dom=this;
                $(dom).DataTable({
			"order": [1],
                	"lengthMenu": [[10, 25, 50, 100, -1], [10, 25, 50, 100, "All"]],
			"pageLength": 10,
		 	"scrollX": "100%",
		});
	//	$(this).parent().toggle();
        });
	$( ".edge-sgout-expand-none" ).parent().next().find('.li-report-content').toggle(false);

	$('.edge-output-direct-datatable').removeClass("ui-table-reflow ui-table");
	$('.edge-output-direct-datatable').each(function(){
		var dom=this;
		var json_table_file = $(this).attr('data-src');
		var dataTableDom = 'lfrtip';
		$.ajax( {
			cache:false,
			url: json_table_file,
			dataType:"json",
			beforeSend: function(){
				//$('#edge-output-datatable-spinner').addClass("edge-sp edge-sp-circle");;
			},
			complete: function(){
			},
			success:function(json){
				var datatable=$(dom).DataTable( {
					"data":json.data,
					"columns":json.columns,
					"dom":dataTableDom,
					"processing": true,
					"deferRender": true,
					"scrollX": true,
					"destroy": true,
					"rowCallback": function( nRow, aData, iDisplayIndex ) {
						if ( aData.Determination== "Positive" ) {
							$(nRow).addClass('edge-fg-red');
						}
						//$('td', nRow).attr('nowrap','nowrap');
						//return nRow;
					},
					"drawCallback": function(settings) {
						$(dom).css({"width":"100%"});
						//contigBlast(tablelinkdom);
					},
					"initComplete": function(settings, json) {
					//	$('#edge-output-datatable-spinner').removeClass("edge-sp edge-sp-circle");;
					}
				});
			},
			error:function(x,t,m){
				//$('#edge-output-datatables-dialog').popup('close');
				showMSG("ACTION FAILED: "+ t + "  Please try again or contact your system administrator." );
			}
		});
		
	});

	var gap_depth_cut_off=0;
	$('.edge-output-datatables').on('click',function(){
		var json_table_file = $(this).attr('data-src');
		var tablelinkdom = $(this);
		var dataTableDom = (/Gap/.test(json_table_file))? '<"#edge-gap-depth-cutoff-div">frtip':'lfrtip';
		$.ajax( {
			cache:false,
			url: json_table_file,
			dataType:"json",
			beforeSend: function(){
				$('#edge-output-datatable-spinner').addClass("edge-sp edge-sp-circle");;
			},
			complete: function(){
			},
			success:function(json){
				var datatable=$('#edge-output-datatable').DataTable( {
					"data":json.data,
					"columns":json.columns,
					"dom":dataTableDom,
					"processing": true,
					"deferRender": true,
					"scrollX": true,
					"destroy": true,
					"rowCallback": function( nRow, aData, iDisplayIndex ) {
						$('td', nRow).attr('nowrap','nowrap');
						return nRow;
					},
					"drawCallback": function(settings) {
						contigBlast(tablelinkdom);
					},
					"initComplete": function(settings, json) {
						$('#edge-output-datatable-spinner').removeClass("edge-sp edge-sp-circle");;
					}
				});
				var gap_depth_input = '<label> Define Gap Depth Coverage: <input type="number" id="edge-gap-depth-cutoff" min="0" max="20" value="'+gap_depth_cut_off+'"> <input type="button" id="edge-gap-depth-cutoff-submit" value="Go"></label>';
				if (/Gap/.test(json_table_file)){
					$('#edge-gap-depth-cutoff-div').html(gap_depth_input);
					$('#edge-gap-depth-cutoff-div').css('display','inline-block');
					$('#edge-gap-depth-cutoff-submit').on('click', function(){
						$('#edge-output-datatables-dialog').popup('close');
						gap_depth_cut_off = $('#edge-gap-depth-cutoff').val();
						defineGapDepth(gap_depth_cut_off,tablelinkdom);
					});
				}
				$('#edge-output-datatables-dialog').popup("reposition",{positionTo: 'window'});
    				$('#edge-output-datatables-dialog').on("popupafterclose", function( event, ui ) {
							datatable.destroy();
							$('#edge-output-datatable').empty();
				});
			},
			error:function(x,t,m){
				$('#edge-output-datatables-dialog').popup('close');
				showMSG("ACTION FAILED: "+ t + "  Please try again or contact your system administrator." );
			}
		});
	});
	function defineGapDepth (cutoff,tablelinkdom){
		$.ajax({
			url: "./cgi-bin/edge_action.cgi",
			dataType: "json",
			cache: false,
			data: { "proj" : focusProjName, "action": 'define-gap-depth', "userType":localStorage.userType, "gap-depth-cutoff" :cutoff,'protocol': location.protocol, 'sid':localStorage.sid},
			beforeSend: function(){
				$.mobile.loading( "show", {
					text: "Doing Gap analysis...",
					textVisible: 1,
					html: ""
				});
			},
			success: function(data){
				$.mobile.loading( "hide" );
				if( data.STATUS == "SUCCESS" ){
					tablelinkdom.attr('data-src',data.PATH);
					tablelinkdom.click();
				}
			},
			error: function(x,t,m){
				$.mobile.loading( "hide" );
				showMSG("ACTION FAILED: "+ t + "  Please try again or contact your system administrator." );
			}
		});
	}

	function contigBlast (tablelinkdom) {
		$('.edge-contigBlast').on('click',function(){
			var contigID=$(this).closest('td').prev('td').text();
			var actionContent = "Do you want to blast " + contigID + " against NT database?<br/>";
			actionContent += "<input type='text' data-mini='true' id='edge-contig-blast-params' placeholder='default: -num_alignments 10 -num_descriptions 10 -evalue 1e-10'>";
			$("#edge_confirm_dialog_content").html(actionContent);
			$('#edge-output-datatables-dialog').popup('close');
			setTimeout(function(){ $( "#edge_confirm_dialog" ).enhanceWithin().popup('open').css('width','360px');},100);
			$("#edge_confirm_dialog a:contains('Cancel')").on("click",function(){
				setTimeout(function(){  tablelinkdom.click();},100);
			});
			$( "#edge_confirm_dialog" ).on('popupafterclose',function(){
				$("#edge_confirm_dialog a").unbind('click');
			});
			$("#edge_confirm_dialog a:contains('Confirm')").unbind('click').on("click",function(){
				var w = window.open("","new");
				w.document.body.innerHTML = '';
				w.document.write( newWindowHeader + "Running Blast on " + contigID + " against NT. Please wait..." + newWindowFooter);
				$.ajax({
					url: "./cgi-bin/edge_action.cgi",
					type: "POST",
					dataType: "json",
					cache: false,
					data: { "proj" : focusProjName, "action": 'contigblast', "userType":localStorage.userType, "contigID" :contigID,'protocol': location.protocol, 'sid':localStorage.sid},
					beforeSend: function(){
						$.mobile.loading( "show", {
							text: "Executing Blast command...",
							textVisible: 1,
							html: ""
						});
					},
					success: function(data){
						data.w = w;
						data.type= 'blast';
						if( data.STATUS == "SUCCESS" ){
							if ( data.PID ){ 
								checkpidInterval = setInterval(function(){check_process(data)},3000); 
							}else{
								$.mobile.loading( "hide" );
								showMSG(data.INFO);
								w.location = edge_path + data.PATH;
								setTimeout(function(){
									w.document.title = newTitle;
									var newHeader= "<div style='background:#50a253;'><h2 style='position:inherit; padding-left:20px;'>" + newTitle + "</h2></div>";
									$(w.document.body).prepend(newHeader);
								},300);
							}
							// add monitor blast file output instead . It may time out for long blast run
						}else{
							$.mobile.loading( "hide" );
							showMSG(data.INFO);
						}
					},
					error: function(x, t, m){
						$.mobile.loading( "hide" );
						setTimeout(function(){ w.close(); },100);
						showMSG("ACTION FAILED: "+ t + "  Please try again or contact your system administrator." );
					},
				});
			});
		});
	}

	$(".edge-piret-report-gene-fc").hide();
	$(":radio[name=edge-piret-output-fc-choose]").on("change", function(){
		if ($(this).val() == "CDS"){
			$(".edge-piret-report-cds-fc").show();
			$(".edge-piret-report-gene-fc").hide();
		}
		if ($(this).val() == "genes"){
			$(".edge-piret-report-cds-fc").hide();
			$(".edge-piret-report-gene-fc").show();
        }
	});
	$(".edge-piret-report-gene-de").hide();
	$(":radio[name^=edge-piret-output-de-choose]").on("change", function(){
		if ($(this).val() == "CDS"){
			$(this).closest(".li-report-content").find(".edge-piret-report-cds-de").show();
			$(this).closest(".li-report-content").find(".edge-piret-report-gene-de").hide();
		}
		if ($(this).val() == "genes"){
			$(this).closest(".li-report-content").find(".edge-piret-report-cds-de").hide();
			$(this).closest(".li-report-content").find(".edge-piret-report-gene-de").show();
        }
	});
	$(".edge-piret-report-fcs-table").removeClass("ui-table-reflow ui-table");
	$(".edge-piret-report-fcs-table").find(".ui-table-cell-label").remove();
	$(".edge-piret-report-fcs-table").DataTable({
		"order": [],
		"deferRender": true,
		"responsive": true,
		"scrollX": "100%",
        });

	$('#edge-output-ccp-rank-by-length-table').hide();
	$(":radio[name=edge-output-ccp-choose-rank]").on("change", function(){
		if($('#edge-output-ccp-choose-rank-length').is(':checked')){
			$('#edge-output-ccp-rank-by-count-table').hide();
			$('#edge-output-ccp-rank-by-length-table').show();
		}
		if($('#edge-output-ccp-choose-rank-count').is(':checked')){
			$('#edge-output-ccp-rank-by-count-table').show();
			$('#edge-output-ccp-rank-by-length-table').hide();
		}
	});

	$('#edge-output-cp-summary').filterable();

	$('#edge-output-cp-choose input').on('click', function(){
		var rank  = $(this).val();
		var rankS = rank.substring(0,3);
		$('#edge-output-cp-summary-filter').val(rank);
		$('#edge-output-cp-summary').filterable("refresh");
		 $('#edge-radarmap-label').html("Radar map at "+ rank +" level ");

		 $(this).closest('li').find('img,iframe,a').each(function(){
			var src = $(this).attr("data-src");
			var asrc = $(this).attr("alt-src");

			if(src){
				var newImg = $(this).clone();
				src = src.replace(/\.(genus|species|strain)\./, "."+rank+"." );
				src = src.replace(/(gen|spe|str)DB/g, rankS+"DB");
				newImg.attr("data-src", src);
			
				if( !$.isEmptyObject(asrc) ){
					asrc = asrc.replace(/\.(genus|species|strain)\./, "."+rank+"." );
					asrc = asrc.replace(/(gen|spe|str)DB/g, rankS+"DB");
					newImg.attr("alt-src", asrc);
				}

				newImg.appendTo($(this).parent());
				newImg.lazyLoadXT();
				$(this).remove();
			}
			else if( $(this).attr("href") ){
				src = $(this).attr("href");
				$(this).attr("href", src.replace(/\.(genus|species|strain)\./, "."+rank+"." ));
			}
		});
	});    

	$('.edge-compare-output').hide()
	$('.edge-compare-species').show();
	$('#edge-compare-rt-choose input').on('click', function(){
		var rank  = $(this).val();
		var rankS = rank.substring(0,3);
		$('.edge-compare-output').hide();
		$('.edge-compare-'+ rank).show();
	});
	$(".edge-compare-output a[data-ajax='false']").each(function(){
		var url=$(this).attr('data-src');
		$(this).attr('href',url);
	});

	var img = $( "#popup-general-img" );
	var obj = $( "#popup-general-obj" );
	$( "#edge-content-report a[data-rel='popup']" ).on( "mouseover", function(){
		var url = $(this).children("img").attr("src");
		if( $(this).children("img").attr("alt-src") ){
			url = $(this).children("img").attr("alt-src");
			$(obj).attr("data", url);
			$(obj).css("background-color", "white");
		}
		else{
			$(img).attr("src", url);
			$(obj).attr("data", url);
		}
	});

	$("table td").tooltipster({multiple:true});

	//krona iframe 
/*	$( "#edge-content-report iframe" ).each(function(){
		var kf = this;		
		var testKronaAnimation = setInterval(function () {
			var h = $("div[title='Help']", $(kf).contents());
			if ( $(h).size() ) {
				clearInterval(testKronaAnimation);
				$(h).parent().hide();
			}
		}, 100);
	});
*/
	//snp tree
	function snptree(url, treetype, svg, name){
		$.get(url, function(data) {
			var dataObject = {
				phyloxml: data,
				fileSource: true
			};
			$("#"+svg).html("");
			var w = $("#"+svg).width() - 10;
			var h = 520;
			
			
			//draw tree
			var name2 = name.replace(/[^a-z0-9_]/gi,'_');
			var phylocanvas = new Smits.PhyloCanvas( dataObject, svg, w, h, treetype );
			$("#"+svg).find("tspan:contains('"+ name2 + "')").css('fill', '#FF0033');
			$("#"+svg).find("tspan:contains('" + name2 + "_contig')").css('fill', 'blue');

			//reset height for rectangular tree
			var otu;
			if( treetype == "rectangular" ){
				otu  = $(dataObject.phyloxml).find("desc").size();
				otu += 2;
				h = 18*otu;
				w = w*1.5;
				var maxX=0;
				var maxY=0;
				var maxTextLen=0;
				var allText = $("#"+svg).find("text");
				for (var idx = 0; idx < allText.length; idx++){
					var textX = parseInt(allText[idx].getAttribute("x"));
					var textY = parseInt(allText[idx].getAttribute("y"));
					if ( maxX < textX){
						maxX = textX;
					}
					if ( maxY < textY){
						maxY = textY;
					}
					if ( maxTextLen < $(allText[idx]).text().length){
						maxTextLen = $(allText[idx]).text().length;
					}
				}
				
				maxX = (maxX * 2 ) + (maxTextLen) * 5;
				maxY = maxY+100;
				$("#"+svg+" > svg" )[0].setAttribute('height', maxY);
				$("#"+svg+" > svg" )[0].setAttribute('width', maxX); 
			}

			//add tooltip
			$("#"+svg+" a").tooltipster();
			$("#"+svg+" a").attr( "target", "new_window");

			if( treetype == "circular" ){
				var svgw  = $("#"+svg+" > svg > path:first-child")[0].getBoundingClientRect().width;
				var svgh = $("#"+svg+" > svg > path:first-child")[0].getBoundingClientRect().height;
				var xoffset = (w-svgw)/2;
				var yoffset = (h-svgh)/2;
				$("#"+svg+" > svg")[0].setAttribute('viewBox',xoffset+" "+yoffset+" "+svgw+" "+svgh)
				$("#"+svg+" > svg" )[0].setAttribute('height', "100%");
				$("#"+svg+" > svg" )[0].setAttribute('width', "100%"); 
			}
		});
	};

	$( "#edge-content-report .output-phylo" ).each(function(){
		var url = $(this).prev().children().first().attr("href");
		var svgid = $(this).attr("id");
		snptree(url, "circular", svgid, projRealname);
	});


	$("#edge-content-report a.sp-circ-btn").on("click", function(){
		var url = $(this).parent().children().first().attr("href");
		var svgid = $(this).parent().next().attr("id");
		snptree(url, "circular", svgid, projRealname);
	});
	
	$("#edge-content-report a.sp-rect-btn").on("click", function(){
		var url = $(this).parent().children().first().attr("href");
		var svgid = $(this).parent().next().attr("id");
		snptree(url, "rectangular", svgid, projRealname);
	});

	$("#edge-content-report a.sp-full-btn").on("click", function(){
		var w = window.open();
		var html = $(this).parent().next().html();
		$(w.document.body).html(html);
		//$(w.document.body).find("svg")
		//	.attr("width", "100%")
		//	.attr("height", "100%");
	});

	//var focusProjName = projdir.split('/').pop();
	// download button
	$("#get_download_link").on('click',function tarProj(){
		var w = window.open("","new","width=360,height=200");
		w.document.body.innerHTML = '';
 		w.document.write(newWindowHeader + "Compressing Project "+ projRealname + " Result ..." + newWindowFooter);
		$.ajax({
			url: "./cgi-bin/edge_action.cgi",
			type: "POST",
			dataType: "json",
			cache: false,
			data: { "proj" : focusProjName, "action": 'tarproj', "userType":localStorage.userType,'protocol': location.protocol, 'sid':localStorage.sid},
			beforeSend: function(){
				$('#get_download_link').addClass('ui-disabled');
				$('#edge-download-spinner').addClass("edge-sp edge-sp-circle");
			},
			complete: function() {
			},
			success: function(data){
				// need data.PATH
				if( data.STATUS == "SUCCESS" ){	
					data.w = w;
					data.spinner_id = "edge-download-spinner";
					data.type = "download";
					if ( data.PID ){ 
						checkpidInterval = setInterval(function(){check_process(data)},3000); 
					}else{
						$('#edge-download-spinner').removeClass("edge-sp edge-sp-circle");
						$('#get_download_link').after(data.LINK);
						$('#ddownload_link').addClass("ui-btn ui-mini ui-btn-inline ui-btn-active");
						$('#ddownload_link').text('Download Project');
						$('#get_download_link').hide();
						showMSG(data.INFO + data.LINK + '.');
						setTimeout(function(){ w.close(); },100);
					}
				}
				else{
					$('#edge-download-spinner').removeClass("edge-sp edge-sp-circle");
					w.close();
					showMSG(data.INFO);
				}
			},
			error: function(data){
				setTimeout(function(){ w.close(); },100);
				$('#edge-download-spinner').removeClass("edge-sp edge-sp-circle");
				showMSG("ACTION FAILED: Please try again or contact your system administrator.");
			}
		});
	});

	var table_row_limit = 10;
	if ($('.edge-get-r2g-reads').length > table_row_limit ){
		$('.edge-get-r2g-reads').closest('table').DataTable({
			"order": [],
			"pageLength": table_row_limit,
			"deferRender": true,
			"responsive": true,
		});
	}
	if ($('.edge-get-c2g-contigs').length > table_row_limit){
		$('.edge-get-c2g-contigs').closest('table').DataTable({
			"order": [],
			"pageLength": table_row_limit,
			"deferRender": true,
			"responsive": true,
		});
	}
	$('.edge-get-c2g-contigs, .edge-get-r2g-reads').on('click',function(){
		var ReferenceID = $(this).closest('tr').find('td:eq(0)').text();
		var ReferenceFile = $(this).closest('tr').find('td:eq(0)').attr("data-reffile");
		var type, action;
		if ($(this).hasClass("edge-get-c2g-contigs")){
			type='contigs';
			action='getcontigbyref';
		}
		else if ($(this).hasClass("edge-get-r2g-reads")){
			type='reads';
			action='getreadsbyref';
		}
		var actionContent = "Do you want to extract mapped to " + ReferenceID + " " + type + " ?<br/>";
		$("#edge_confirm_dialog_content").html(actionContent);
		$( "#edge_confirm_dialog" ).enhanceWithin().popup('open').css('width','360px');
		$("#edge_confirm_dialog a:contains('Confirm')").unbind('click').on("click",function(){
			var w = window.open("","new","width=360,height=200");
			w.document.body.innerHTML = '';
			w.document.write( newWindowHeader + "Extracting Contigs/Reads Mapped to " + ReferenceID + ". Please wait..." + newWindowFooter);	
			$.ajax({
				url: "./cgi-bin/edge_action.cgi",
				type: "POST",
				dataType: "json",
				cache: false,
				data: { "proj" : focusProjName, "action": action, "reffile":ReferenceFile,"refID":ReferenceID,"userType":localStorage.userType,'protocol': location.protocol, 'sid':localStorage.sid},
				beforeSend: function(){
					$.mobile.loading( "show", {
						text: "Extract Contigs/Reads Fasta/Fastq...",
						textVisible: 1,
						html: ""
					});
				},
				complete: function() {
				},
				success: function(data){
					if( data.STATUS == "SUCCESS" ){
						data.w = w;
						if ( data.PID ){ 
							checkpidInterval = setInterval(function(){check_process(data)},3000); 
						}else{
							$.mobile.loading( "hide");
                                        		//console.log(edge_path,data.PATH);
							w.opener.location = edge_path + data.PATH;
							setTimeout(function(){ w.close(); },100);
						}
					}else{
						$.mobile.loading( "hide");
						w.close();
						showMSG(data.INFO);
					}
				},
				error: function(data){
					$.mobile.loading( "hide");
					setTimeout(function(){ w.close(); },100);
					showMSG("ACTION FAILED: Please try again or contact your system administrator.");
				}
			});
		});
	});
	$('#edge-get-contigs-by-taxa').on('change',function(){
		var taxa = $(this).val();
		if(taxa == "0" ){
			return;
		}
		var w = window.open("","new","width=360,height=200");
		w.document.body.innerHTML = '';
		w.document.write(newWindowHeader + "Extracting " + taxa + " FASTA. Please wait ..." + newWindowFooter);
		$.ajax({
			url: "./cgi-bin/edge_action.cgi",
			type: "POST",
			dataType: "json",
			cache: false,
			data: { "proj" : focusProjName, "action": 'getcontigbytaxa', "taxa":taxa,"userType":localStorage.userType,'protocol': location.protocol, 'sid':localStorage.sid},
			beforeSend: function(){
				$('#edge-get-contigs-by-taxa').addClass('ui-disabled');
				$('#edge-get-contigs-spinner').addClass("edge-sp edge-sp-circle");
			},
			complete: function() {
			},
			success: function(data){
				if( data.STATUS == "SUCCESS" ){
					data.spinner_id = 'edge-get-contigs-spinner';
					data.w = w;
					data.select_id = 'edge-get-contigs-by-taxa';
					if ( data.PID ){ 
						checkpidInterval = setInterval(function(){check_process(data)},3000); 
					}else{
						$('#edge-get-contigs-by-taxa').removeClass('ui-disabled');
						$('#edge-get-contigs-spinner').removeClass("edge-sp edge-sp-circle");
						//console.log(edge_path,data.PATH);
						w.opener.location = edge_path + data.PATH;
						setTimeout(function(){ w.close(); },100);
					}
				}else{
					$('#edge-get-contigs-spinner').removeClass("edge-sp edge-sp-circle");
					w.close();
					showMSG(data.INFO);
				}
			},
			error: function(data){
				$('#edge-get-contigs-spinner').removeClass("edge-sp edge-sp-circle");
				setTimeout(function(){ w.close(); },100);
				showMSG("ACTION FAILED: Please try again or contact your system administrator.");
			}
		});
		
	});
	$('.edge-get-reads-by-taxa').on('change',function(){
		var select_id = this.id;
		var taxa = $('#'+ select_id + ' option:selected ').val();
		var taxa_name = $('#'+ select_id + ' option:selected ').text();
		var cptool = select_id.replace('edge-get-reads-by-taxa-','');
		var w = window.open("","new","width=360,height=200");
		w.document.body.innerHTML = '';
		w.document.write( newWindowHeader  + "Extracting " + taxa_name + " FASTQ from " + cptool + ". Please wait..." + newWindowFooter);
		$.ajax({
			url: "./cgi-bin/edge_action.cgi",
			type: "POST",
			dataType: "json",
			cache: false,
			data: { "proj" : focusProjName, "action": 'getreadsbytaxa', "cptool":cptool,"taxa":taxa,"userType":localStorage.userType,'protocol': location.protocol, 'sid':localStorage.sid},
			beforeSend: function(){
				$('#'+ select_id).addClass('ui-disabled');
				$('#edge-get-reads-spinner-'+cptool).addClass("edge-sp edge-sp-circle");
				$('#edge-get-reads-spinner-'+cptool).css("width","20px").css("height","20px");
			},
			complete: function(data) {
			},
			success: function(data){
				//console.log(data);
				if( data.STATUS == "SUCCESS" ){
					data.spinner_id = 'edge-get-reads-spinner-'+cptool;
					data.w = w;
					data.select_id = select_id;
					if ( data.PID ){ 
						checkpidInterval = setInterval(function(){check_process(data)},3000); 
					}else{
						$('#'+ select_id).removeClass('ui-disabled');
						$('#edge-get-reads-spinner-'+cptool).removeClass("edge-sp edge-sp-circle");
						w.opener.location = edge_path + data.PATH;
						setTimeout(function(){ w.close(); },300);
					}
				}else{
					$('#edge-get-reads-spinner-'+cptool).removeClass("edge-sp edge-sp-circle");
					w.close();
					showMSG(data.INFO);
				}
			},
			error: function(data){
				$('#edge-get-reads-spinner-'+cptool).removeClass("edge-sp edge-sp-circle");
				setTimeout(function(){ w.close(); },100);
				showMSG("ACTION FAILED: Please try again or contact your system administrator.");
			}
		});
		
	});

	function check_process(data){
		var spinner_id = data.spinner_id;
		var w = data.w;
		var select_id = data.select_id;
		$.ajax({
			url: "./cgi-bin/edge_action.cgi",
			type: "POST",
			dataType: "json",
			cache: false,
			data: { "action": 'checkpid', "pid": data.PID,'protocol': location.protocol, 'sid':localStorage.sid},
			success: function(obj){
				if( obj.STATUS == "DONE" ){
					clearInterval(checkpidInterval);
					if (data.type === "download"){
 						$('#get_download_link').after(data.LINK);
 						$('#ddownload_link').addClass("ui-btn ui-mini ui-btn-inline ui-btn-active");
 						$('#ddownload_link').text('Download Project');
						$('#get_download_link').hide();
 						showMSG(data.INFO + data.LINK + '.');
					}
					if (spinner_id){
						$('#'+ select_id).removeClass('ui-disabled');
						$('#' + spinner_id).removeClass("edge-sp edge-sp-circle");
					}else{
						$.mobile.loading( "hide");
					}
					
					if (data.type === "blast"){
						showMSG(data.INFO);
						w.location = edge_path + data.PATH;
						setTimeout(function(){
							w.document.title = newTitle;
							var newHeader= "<div style='background:#50a253;'><h2 style='position:inherit; padding-left:20px;' class='edge-header'>"+newTitle+"</h2></div>";
							$(w.document.body).prepend(newHeader);
						},300);
					}else{
						w.opener.location = edge_path + data.PATH;
						setTimeout(function(){ w.close(); },300);
					}
				}else{
					//w.document.write(".");
					//console.log(obj.INFO);
				}
			},
			error: function(obj){
				showMSG("ACTION FAILED: Please try again or contact your system administrator.");
			}
		});
		
	}
	function showMSG( dialog_content ) {
		$( "#edge_integrity_dialog_header" ).text("Message");
		$( "#edge_integrity_dialog_content" ).html(dialog_content);
		$( "#edge_integrity_dialog" ).popup('open');
	}
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
	// Take Notes actions
	var NotesLenLimit = 500,
		NotesLen; // Maximum word length
	$('#edge-takenotes-content').keydown(function(event) {	
		NotesLen = $('#edge-takenotes-content').val().split(/[\s]+/);
		if (NotesLen.length >= NotesLenLimit) { 
			if ( event.keyCode == 46 || event.keyCode == 8 ) {// Allow backspace and delete buttons
			} else if (event.keyCode < 48 || event.keyCode > 57 ) {//all other buttons
    				event.preventDefault();
    			}
		}
		//console.log(NotesLen.length + " words are typed out of an available " + NotesLenLimit);
		var wordsLeft = (NotesLenLimit) - NotesLen.length;
		$('.NotesWords-left').html(wordsLeft + ' words left');
	});
	$('#edge-result-notes').on('click',function(){
		$( "#edge-takenotes-dialog" ).popup('open');
	});
	$( "#edge-takenotes-dialog" ).popup({
		beforeposition: function( event, ui ) {
                        getNotes(projNotes);
                }
	});
	function getNotes(texturl) {
		$.ajax({
			url: texturl,
			dataType: 'text',
			cache: false,
			success: function(text) {
				$("#edge-takenotes-content").val(text);
			},
			error: function(text) {
				$("#edge-takenotes-content").attr("placeholder","Enter Note Here...(alphanumeric space comma dot dash underscore) "); 
			}
		});
	}
	$("#edge-takenotes-dialog a:contains('Save')").on("click",function(){
		var text = $("#edge-takenotes-content").val();
		$.ajax({
			url: "./cgi-bin/edge_action.cgi",
			type: "POST",
			dataType: "json",
			cache: false,
			data: { "action": 'takenotes', "projnotes": text, "proj" : focusProjName, 'protocol': location.protocol, 'sid':localStorage.sid},
			success: function(data){
				if (data.INFO){
					showMSG(data.INFO);
				}
			},
			error: function(data){
				showMSG(data.INFO);
				//showMSG("ACTION FAILED: Please try again or contact your system administrator.");
			}
		});
	});
        $('.edge-output-tool-summary').removeClass("ui-table-reflow ui-table");
	$('.edge-output-tool-summary').find(".ui-table-cell-label").remove();
	$('.edge-output-tool-summary').each(function(){
		var dom=this;
		$(dom).DataTable({
			"pageLength": 10,
		});
	});
	$(".rtaxToolresult").hide();
	var firstTaxToolResultTag = $('#showTaxSelect option').eq(1).val();
	$('#'+ firstTaxToolResultTag).show();
	$('#showTaxSelect').on('change',function(){
		$(".rtaxToolresult").hide();
		$('#showTaxSelect option:selected').each(function(){
			var tag=$(this).val();
			$('#'+ tag).show();
		});
	});

	$('.showTax').first().click();
});
//# sourceURL=edge-output.js
