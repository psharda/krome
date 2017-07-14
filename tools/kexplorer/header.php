<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<title>KROME rate explorer</title>
<link href='http://fonts.googleapis.com/css?family=Open+Sans' rel='stylesheet' type='text/css' />
<link href="css/styles.css" rel="stylesheet" type="text/css" />
<link rel="stylesheet" href="jquery-ui-1.12.1.custom/jquery-ui.css">
<link rel="stylesheet" href="jquery-ui-1.12.1.custom/style.css">
<script src="jquery-ui-1.12.1.custom/external/jquery/jquery.js"></script>
<script src="jquery-ui-1.12.1.custom/jquery-ui.js"></script>
<script>

	function getColor(value){
	    //value from 0 to 1
	    var hue=((1-value)*120).toString(10);
	    return ["hsl(",hue,",50%,50%)"].join("");
	}

	//load JSON
	function loadJSON(idx){
		//parse JSON
		$.getJSON("data.json", function(json) {
			hline = "<tr><th><th><th><th><th><th><th><th><th><th><th><th><th><th><th>"
			//create a table with fluxes
			table = "<table width=\"100%\">"
			table += hline
			table += "<tr><td>#<td><td><td><td><td><td><td><td><td><td><td><td><td>flux cm<sup>-3</sup>s<sup>-1</sup><td>fTot*"
			table += hline
			for (i=0; i<json[idx]["data"].length; i++){
				pFlux = json[idx]["data"][i]["fluxNormTot"]
				cFlux = 1.-Math.max(0e0,Math.log10(pFlux+1e-40)+8.)/8.
				table += "<tr style=\"background-color:"+getColor(cFlux)+";\"><td>"+(i+1)+".<td>"+json[idx]["data"][i]["htmlVerbatim"]
					+ "<td>"+(json[idx]["data"][i]["flux"]).toExponential(2)
					+ "<td>"+pFlux.toExponential(5)
			}
			table += hline
			table += "</table><br>"
			table += "*fraction of total flux"

			//populate html spans
			$("#slider-value").html(table);
			$("#slider-img").html("<img src=\"pngs/graph_"+(1e6+idx)+".png\">");
			$("#slider-xvar").html(json[idx]["xvar"].toExponential(4));
		});
	}

	//init JQuery slider
	$(function() {
		$("#slider").slider(
			{
			    value:0,
			    min: 0,
			    max: #MAXFILE#,
			    step: 1,
			    slide: function( event, ui ) {
				loadJSON(ui.value);
			    }
			}
		);
	});
</script>

</head>


<body onload="loadJSON(0);">
<div class="header-wrapper">
  <div class="header">
    <div class="logo">
      <h1>KROME explorer</h1>
    </div>
    <!---logo ends here -->
    <!--<div class="menu">
      <ul>
        <li><a class="active" href="index.html">Home</a>me</li>
      </ul>
    </div>-->
    <!--- menu ends here -->
  </div>
  <!---headerends here -->
</div>

<div class="clear"></div>
<!---container ends here-->
<div class="conainer">


