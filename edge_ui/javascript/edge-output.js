$( document ).ready(function()
{
	$( "#edge-out-expand-all" ).on("click", function(){
		$( "#edge-content-report > div[data-role='collapsible']" ).collapsible( "option", "collapsed", false );
	})

	$( "#edge-out-expand-none" ).on("click", function(){
		$( "#edge-content-report > div[data-role='collapsible']" ).collapsible( "option", "collapsed", true );
	})

	$( "#edge-content-report p.extra-text-link > a" ).attr( "target", "new_window");
	$( "#edge-content-report .iframe_label > a" ).attr( "target", "new_window");
	$( "#edge-content-report .jb_link" ).attr( "target", "new_window");
	$( "#edge-content-report table.output-table td > a" ).attr( "target", "new_window");

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
	$('.edge-outputfile-tree').removeAttr("target");
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
		}
	});

	$("table td").tooltipster();

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
			
			//reset height for rectangular tree
			var otu;
			if( treetype == "rectangular" ){
				otu  = $(dataObject.phyloxml).find("desc").size();
				otu += 2;
				h = 14*otu;
				w = w*1.5;
			}
			
			//draw tree
			var phylocanvas = new Smits.PhyloCanvas( dataObject, svg, w, h, treetype );
			$("#"+svg).find('tspan:contains("_contig")').css('fill', 'blue');
			$("#"+svg).find('tspan:contains("_read")').css('fill', 'orange');

			//add tooltip
			$("#"+svg+" a").tooltipster();
			$("#"+svg+" a").attr( "target", "new_window");

			if( treetype == "circular" ){
				var svgw  = $("#"+svg+" > svg > path:first-child")[0].getBoundingClientRect().width;
				var svgh = $("#"+svg+" > svg > path:first-child")[0].getBoundingClientRect().height;
				var xoffset = (w-svgw)/2;
				var yoffset = (h-svgh)/2;
				$("#"+svg+" > svg")[0].setAttribute('viewBox',xoffset+" "+yoffset+" "+svgw+" "+svgh)
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
		$(w.document.body).find("svg")
			.attr("width", "100%")
			.attr("height", "100%");
	});

});

//# sourceURL=edge-output.js
