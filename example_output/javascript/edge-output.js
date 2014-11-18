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

	var img = $( "#popup-general-img" );
	$( "#edge-content-report a[data-rel='popup']" ).on( "mouseover", function(){
		var url = $(this).children("img").attr("src");
		$(img).attr("src", url);
	});

	$("table td").tooltipster();

	//krona iframe 
	$( "#edge-content-report iframe" ).each(function(){
		var kf = this;		
		var testKronaAnimation = setInterval(function () {
			var h = $("div[title='Help']", $(kf).contents());
			if ( $(h).size() ) {
				clearInterval(testKronaAnimation);
				$(h).parent().hide();
			}
		}, 100);
	});

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
			phylocanvas = new Smits.PhyloCanvas( dataObject, svg, w, h, treetype );
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

	//display every graphics
	$( "#edge-content-report iframe" ).show();
	$( "div[class^='ui-grid-']" ).show();
});

//# sourceURL=edge-output.js
