$( document ).ready(function()
{
	var projdir = $('#edge-proj-outdir').attr('dir-src');
	var focusProjName = $('#edge-output-projname').attr("data-pid");

	$( "#edge-out-expand-all" ).on("click", function(){
		$( "#edge-content-report > div[data-role='collapsible']" ).collapsible( "option", "collapsed", false );
	});

	$( "#edge-out-expand-none" ).on("click", function(){
		$( "#edge-content-report > div[data-role='collapsible']" ).collapsible( "option", "collapsed", true );
	});

	$( "#edge-content-report p.extra-text-link > a:not(a[data-rel='popup'])" ).attr( "target", "new_window");
	$( "#edge-content-report .iframe_label > a" ).attr( "target", "new_window");
	$( "#edge-content-report .jb_link" ).attr( "target", "new_window");
	$( "#edge-content-report table.output-table td > a:not(a[data-rel='popup'])" ).attr( "target", "new_window");
	$( ".edge-compare-output .iframe_label > a" ).attr( "target", "_blank");

	$( "#edge-content-report span.li-report-content-title" ).on("click",function(){
		$(this).next().slideToggle();
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

	// directory file tree
	var dir;
	var loc = window.location.pathname;
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
	var json_table_file;
	$('.edge-output-datatables').on('click',function(){
		json_table_file = $(this).attr('data-src');
		var tablelinkdom = $(this);
		$.ajax( {
			cache:false,
			url: json_table_file,
			dataType:"json",
			beforeSend: function(){
			},
			complete: function(){
			},
			success:function(json){
				var datatable=$('#edge-output-datatable').dataTable( {
					"aaData":json.data,
					"aoColumns":json.columns,
					//"deferRender": true,
					"scrollX": true,
					"destroy": true,
					"bDestroy":true,
					"fnRowCallback": function( nRow, aData, iDisplayIndex ) {
						$('td', nRow).attr('nowrap','nowrap');
						return nRow;
					},
					"drawCallback": function(settings) {
						contigBlast(tablelinkdom);
					}
				});
				$('#edge-output-datatables-dialog').popup("reposition",{positionTo: 'window'});
    			$('#edge-output-datatables-dialog').on("popupafterclose", function( event, ui ) {
							datatable.fnDestroy();
							$('#edge-output-datatable').empty();
				});
			}
		});
	});

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
						$.mobile.loading( "hide" );
						if( data.STATUS == "SUCCESS" ){
							showMSG(data.INFO);
							window.open(data.PATH);
							// add monitor blast file output instead . It may time out for long blast run
						}else{
							showMSG(data.INFO);
						}
					},
					error: function(x, t, m){
						$.mobile.loading( "hide" );
						showMSG("ACTION FAILED: "+ t + "  Please try again or contact your system administrator." );
					},
				});
			});
		});
	}

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
	$('#edge-compare-summary-table').find('a').each(function(){
		var projId=$(this).attr('data-pid');
		var url = "/?proj=" + projId;
		$(this).attr('href',url);
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
	function snptree(url, treetype, svg){
		$.get(url, function(data) {
			var dataObject = {
				phyloxml: data,
				fileSource: true
			};
			$("#"+svg).html("");
			var w = $("#"+svg).width() - 10;
			var h = 520;
			
			
			//draw tree
			var phylocanvas = new Smits.PhyloCanvas( dataObject, svg, w, h, treetype );
			$("#"+svg).find('tspan:contains("_contig")').css('fill', 'blue');
			$("#"+svg).find('tspan:contains("read")').css('fill', '#FF0033');

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
		snptree(url, "circular", svgid);
	});


	$("#edge-content-report a.sp-circ-btn").on("click", function(){
		var url = $(this).parent().children().first().attr("href");
		var svgid = $(this).parent().next().attr("id");
		snptree(url, "circular", svgid);
	});
	
	$("#edge-content-report a.sp-rect-btn").on("click", function(){
		var url = $(this).parent().children().first().attr("href");
		var svgid = $(this).parent().next().attr("id");
		snptree(url, "rectangular", svgid);
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
				if( data.STATUS == "SUCCESS" ){
					$('#edge-download-spinner').removeClass("edge-sp edge-sp-circle");
					$('#get_download_link').after(data.LINK);
					$('#ddownload_link').addClass("ui-btn ui-mini ui-btn-inline ui-btn-active");
					$('#ddownload_link').text('Download Project');
					$('#get_download_link').hide();
					showMSG(data.INFO + data.LINK + '.');
				}
				else{
					$('#edge-download-spinner').removeClass("edge-sp edge-sp-circle");
					showMSG(data.INFO);
				}
			},
			error: function(data){
				$('#edge-download-spinner').removeClass("edge-sp edge-sp-circle");
				showMSG("ACTION FAILED: Please try again or contact your system administrator.");
			}
		});
	});
	$('#edge-get-contigs-by-taxa').on('change',function(){
		var taxa = $(this).val();
		if(taxa == "0" ){
			return;
		}
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
					$('#edge-get-contigs-by-taxa').removeClass('ui-disabled');
					$('#edge-get-contigs-spinner').removeClass("edge-sp edge-sp-circle");
					//console.log(edge_path,data.PATH);
					window.open(edge_path + data.PATH);
				}else{
					$('#edge-get-contigs-spinner').removeClass("edge-sp edge-sp-circle");
					showMSG(data.INFO);
				}
			},
			error: function(data){
				$('#edge-get-contigs-spinner').removeClass("edge-sp edge-sp-circle");
				showMSG("ACTION FAILED: Please try again or contact your system administrator.");
			}
		});
		
	});
	$('.edge-get-reads-by-taxa').on('change',function(){
		var select_id = this.id;
		var taxa = $('#'+ select_id).val();
		var cptool = select_id.replace('edge-get-reads-by-taxa-','');
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
			complete: function() {
			},
			success: function(data){
				if( data.STATUS == "SUCCESS" ){
					$('#'+ select_id).removeClass('ui-disabled');
					$('#edge-get-reads-spinner-'+cptool).removeClass("edge-sp edge-sp-circle");
					window.open(edge_path + data.PATH);
				}else{
					$('#edge-get-reads-spinner-'+cptool).removeClass("edge-sp edge-sp-circle");
					showMSG(data.INFO);
				}
			},
			error: function(data){
				$('#edge-get-reads-spinner-'+cptool).removeClass("edge-sp edge-sp-circle");
				showMSG("ACTION FAILED: Please try again or contact your system administrator.");
			}
		});
		
	});
	function showMSG( dialog_content ) {
		$( "#edge_integrity_dialog_header" ).text("Message");
		$( "#edge_integrity_dialog_content" ).html(dialog_content);
		$( "#edge_integrity_dialog" ).popup('open');
	}
});
//# sourceURL=edge-output.js
