function addHostList(){
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
}

var target_menu = "#edge-ref-file-fromlist-menu,#edge-phylo-ref-select-menu,#edge-hostrm-file-fromlist-menu,#edge-get-contigs-by-taxa-meun,.edge-get-reads-by-taxa";
var target_dialog = ["edge-ref-file-fromlist-dialog","edge-phylo-ref-select-dialog","edge-hostrm-file-fromlist-dialog","edge-get-contigs-by-taxa-dialog"];

$.mobile.document
    // The custom selectmenu plugin generates an ID for the listview by suffixing the ID of the
    // native widget with "-menu". Upon creation of the listview widget we want to place an
    // input field before the list to be used for a filter.
    .on( "listviewcreate", target_menu, function( event ) {
        var input,
            list = $( event.target ),
            form = list.jqmData( "filter-form" ),
            id = this.id.replace("-menu",""),
            $select_menu = $( "#" + id);
        // We store the generated form in a variable attached to the popup so we avoid creating a
        // second form/input field when the listview is destroyed/rebuilt during a refresh.
        if ( !form ) {
            input = $( "<input data-type='search' placeholder='ex: Escherichia'></input>" );
			select = $("<div id='" + id + "-allnone-btn' style='width:130px;margin-left:auto;padding:0.3em'> Select <a href='#' id='"+ id + "-all'>All</a> | <a href='#' id='" + id + "-none'>None</a> </div>");
			form = $( "<form></form>" ).append( input ).append(select);
            //form = $( "<form></form>" ).append( input );
            input.textinput();
            list
                .before( form )
                .jqmData( "filter-form", form ) ;
            form.jqmData( "listview", list );
        }
        // Instantiate a filterable widget on the newly created listview and indicate that the
        // generated input form element is to be used for the filtering.
        if (id.indexOf('hostrm')>0 || id.indexOf('get-contigs')>0 || id.indexOf('get-reads')>0 ){
        	list.filterable({
        		input: input,
            	children: "> li:not(:jqmData(placeholder='true'))",
            	filterCallback: OrSearch
        	});
        }else{
			list.filterable({
				input: input,
				children: "> li:not(:jqmData(placeholder='true'))",
				filterCallback: function (){},
				beforefilter: function ( e, data ) {
					var $ul = $( this ),
					$input = $( data.input ),
					value = $input.val(),
					idx =0;
					$ul.html( "" );
					$select_menu.html("");
					var queries = value.split(/[ ,]+/);
					var query="";
					for ( idx = 0 ; idx < queries.length ; idx++ ) {
							if (queries[idx] && queries[idx].length > 3){
								query += queries[idx] + " ";
							}
					}
				  
					if ( query && query.length > 3 ) {
						$ul.html( "<li><div class='ui-loader'><span class='ui-icon ui-icon-loading'></span></div></li>" );
						$ul.listview( "refresh" );
						$.ajax({
						  url: "./cgi-bin/edge_dynamic_selection.cgi",
						  dataType: "json",
						  data: { query: query }
						  })
						  .then( function ( response ) {
								var genomeIds = Object.keys(response);
								genomeIds.sort();
								var testToInsert = [];
								var OptiontToInsert = [];
								var index = 0;
								while ( index < genomeIds.length){
									var name= genomeIds[index];
									name = name.replace(/_uid\d+$/, "");
									OptiontToInsert[index] = "<option value="+name+">"+name+"</option>";
									testToInsert[index] = "<li>" + name + "</li>";
									index++;
								}
								$ul.html( testToInsert.join('') );
								$ul.listview( "refresh" );
								//$ul.trigger( "updatelayout");
								$select_menu.html( "<option disabled='' data-placeholder='true'>Select reference genome(s)...</option>" + OptiontToInsert.join('') );
								$select_menu.selectmenu( "refresh" );
								//$select_menu.trigger( "updatelayout");
							
							
						});
					}
				},
				filter: function( event, ui ) {
						if($(this).children().length){
							$(this).children().first().addClass("ui-screen-hidden");
							$(this).find('a').each(function(){
								var wordLen = $(this).html().length;
								if (wordLen> 50){
									$(this).css("font-size","11px");
								}
								else if (wordLen> 45){
									$(this).css("font-size","12px");
								}
								else if (wordLen> 40){
									$(this).css("font-size","13px");
								}else {
									$(this).css("font-size","14px");
								}
							});
						//	console.log($(this).children().first());
						}else{
							//$(this).html( "<li>Select reference genome(s)...</li>" );
						//	$(this).listview( "refresh" );
						//	$select_menu.html( "<option disabled='' data-placeholder='true'>Select reference genome(s)...</option>");
						//	$select_menu.selectmenu( "refresh" );
						}
				}
			 });
         }
        // select all or none 
         $("#"+id+"-none").click(function() {
   	     	list.children().not(".ui-screen-hidden").children("a").removeClass("ui-checkbox-on");
   	     	list.children().not(".ui-screen-hidden").children("a").addClass("ui-checkbox-off");
   	     	$select_menu.find("option:not([disabled])").removeAttr("selected");
   	     	$select_menu.selectmenu('refresh');
   	 });
   	 $("#"+id+"-all").click(function() {
   	     	list.children().not(".ui-screen-hidden").children("a").removeClass("ui-checkbox-off");	
   	     	list.children().not(".ui-screen-hidden").children("a").addClass("ui-checkbox-on");
   	     	$select_menu.find("option:not([disabled])").attr("selected",'selected');
		list.children().not(".ui-screen-hidden").attr("aria-selected",true);
   	     	$select_menu.selectmenu('refresh');
   	 });
    })
    // The custom select list may show up as either a popup or a dialog, depending on how much
    // vertical room there is on the screen. If it shows up as a dialog, then the form containing
    // the filter input field must be transferred to the dialog so that the user can continue to
    // use it for filtering list items.
    .on( "pagecontainerbeforeshow", function( event, data ) {
        var listview, form,
            id = data.toPage && data.toPage.attr( "id" );
        if ( !( id && (target_dialog.indexOf(id)>=0 || id.indexOf('get-reads')>0) ))
	{
            return;
        }
        listview = data.toPage.find( "ul" );
        form = listview.jqmData( "filter-form" );
        // Attach a reference to the listview as a data item to the dialog, because during the
        // pagecontainerhide handler below the selectmenu widget will already have returned the
        // listview to the popup, so we won't be able to find it inside the dialog with a selector.
        data.toPage.jqmData( "listview", listview );
        // Place the form before the listview in the dialog.
        listview.before( form );
    })
    // After the dialog is closed, the form containing the filter input is returned to the popup.
    .on( "pagecontainerhide", function( event, data ) {
        var listview, form,
            id = data.toPage && data.toPage.attr( "id" );
        if ( !( id && (target_dialog.indexOf(id)>=0 || id.indexOf('get-reads')>0) ))
	{
            return;
        }
        listview = data.toPage.jqmData( "listview" );
	if ( ! listview ) {
        	listview = data.prevPage.jqmData( "listview" );
        }
        popupID = id.replace("-dialog","-listbox-popup");
        form = listview.jqmData( "filter-form" );
	if ( !form ) {
		var input = $( "<input data-type='search'></input>" );
		form = $( "<form></form>" ).append( input );
		listview.filterable({
		  input: input,
		  children: "> li:not(:jqmData(placeholder='true'))",
		  filterCallback: OrSearch
		});
		listview.jqmData( "filter-form",form);
	}
        // Put the form back in the popup. It goes ahead of the listview.
        listview.before( form );
	$('#'+popupID).on("popupbeforeposition", function( event, ui ) {
        	listview.before( form );
	})
	AddSelectRefList();
	$('#'+id).enhanceWithin();
    });

$( document ).ready(function(){
	//addHostList();
	$('#edge-phylo-ref-select-listbox').on( "popupafterclose", function(){
		AddSelectRefList();
	});
	$('#edge-hostrm-file-fromlist-listbox').on( "popupafterclose", function(){
		var $selected = $("#edge-hostrm-file-fromlist option:selected");
		if ($selected.size() > 0){
			$('#edge-hostrm-sw1').click().checkboxradio("refresh");
		}
	});
	
});

function AddSelectRefList() {
	var selectedList = "<option value=0 selected='selected'>Random</option>";
	$selected = $("#edge-phylo-ref-select option:selected");
	if ($selected.size() > 0){
		$("#edge-phylo-ref-select-ref-div").show();  	
		$selected.each( function( key, val ) {
			var value = $(this).html();
			selectedList += "<option value=" + value + ">" + value + "</option>"; 
		});
		$("#edge-phylo-ref-select-ref").html(selectedList);
		$("#edge-phylo-ref-select-ref").selectmenu( "refresh" );
   	 }else{
   	 	$("#edge-phylo-ref-select-ref-div").hide();
   	 }
}   

function OrSearch( index, searchValue ) {
	var ret = false;
	var idx;
	var notMatch_cnt=0;
	var search_cnt=0;
	if (searchValue && searchValue.length > 3){  
		searchValue = searchValue.split(/[ ,]+/);
	    var text = $(this).text().toLowerCase();
	    var filttext = $(this).data("filtertext") || '';
	    filttext = filttext.toLowerCase();
		for ( idx = 0 ; idx < searchValue.length ; idx++ ) {
			if (searchValue[idx] && searchValue[idx].length > 0){
				search_cnt++;
				if( text.indexOf(searchValue[idx].toLowerCase(), 0) < 0 && filttext.indexOf(searchValue[idx].toLowerCase(), 0) < 0){
					notMatch_cnt++;
				}
			}
	    }
	    if (search_cnt === notMatch_cnt) {
			ret = true;
	    }
	} 
	return ret;
}
