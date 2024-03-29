!(function (d3) {

$("bcontent").empty();


    // manhattan plot code obtained/edited/adapted from : https://github.com/chengsoonong/rede/tree/master/Rede/manhattanplot
    var chartDiv = document.getElementById("containergenome2");
    d3.select("#containergenome2");
//                .style("height", "100%")
//     .style("width", "100%");
//                .style("overflow", "auto");

    var chrom_length = 0;

	var chrom_acum_length = new Array();
	var chromLength = new Array(249250621, 243199373, 198022430, 191154276,
    180915260, 171115067, 159138663, 146364022,
    141213431, 135534747, 135006516, 133851895,
    115169878, 107349540, 102531392, 90354753,
    81195210, 78077248, 59128983, 63025520,
    48129895, 51304566, 155270560, 59373566);

	//this initializes chrom_length and chrom_acum_length to be used in manhattan plot
	for (var i = 0; i < chromLength.length; i++) {
  	  chrom_length = chrom_length + chromLength[i];
   	 chrom_acum_length.push(chrom_length);
	}


	posarr = new Array();
	data_low = new Array();
        data_high = new Array();
 	pvalarr = new Array();
        idarr = new Array(); 
	data_weight_pvalue = [];
        var cutoff = 1;
    d3.csv("/static/mrpgene/phebf/site/{{namespace}}/" + icd[0].icd + ".chrpos.01_0001.ptv.out",function(error,data){
    data.forEach(function(d) {
	// Fill two seperate arrays, one with points greater than the -log10(pvalue) threshold, and another
	// with points lower than the threshold defined in cutoff
	       d.chrom = +d.chrom;
	       d.pos = +d.pos;
			    //chrom,pos,gene,l10bfindep,l10bfsim,nvar 
	       d.l10bf = +d.l10bfindep;
         	if (d.chrom == 1) {
                if (d.l10bf < cutoff) {
                        data_low.push([d.pos, d.l10bf,
                    	d.gene, d.chrom]);
                } else {
					data_high.push([d.pos, d.l10bf,
					d.gene, d.chrom]);
				}
			} else {
				if (d.l10bf < cutoff) {
					data_low.push([d.pos + chrom_acum_length[d.chrom - 2],
						d.l10bf,
						d.gene, d.chrom]);
				} else {
						data_high.push([d.pos + chrom_acum_length[d.chrom - 2],
							d.l10bf,
							d.gene, d.chrom]);


	}			

        }
        		        	data_weight_pvalue.push(d.l10bf);
    
             
		})

	// Merge the both high and low significance points together into a single array
    all_data = data_low.concat(data_high);
    var min_pvalue = d3.max([-5, d3.min([-5, d3.min(data_weight_pvalue)])]);
    var max_pvalue = d3.max([5, d3.max(data_weight_pvalue)]);
    // var which defined the extra space in the bottom and the top of the -log(p-value) axis dynamically to the dataset
    var extend_scale = (max_pvalue - min_pvalue) *0.05;
    ix_1 = 0;
    ix_2 = chrom_length;
    iy_1 = min_pvalue - extend_scale;
    iy_2 = max_pvalue + extend_scale;

	//create the manhattan plot
	x1 = ix_1;
	x2 = ix_2;
	y1 = iy_1;
	y2 = iy_2;

    var margin = {
        top: 30,
        right: 50,
        bottom: 20,
        left: 60
    };
    var colorScale = d3.scale.log()
        .domain([d3.min(all_data, function(d) {
            return parseInt(d[3]);
        }), d3.max(all_data, function(d) {
            return parseInt(d[3]);
        })])
        .interpolate(d3.interpolateHsl)
        .range(["#B3995D", "#8C1515"]);

    var w = chartDiv.clientWidth - margin.left - margin.right; //900;
    var h = 500 - margin.top - margin.bottom; //600;
	console.log(w,h,x1,x2,y1,y2);

    //Create scale functions
    var xScale = d3.scale.linear()
        .domain([x1, x2])
        .range([0, w]);

    var yScale = d3.scale.linear()
        .domain([y1, y2])
        .range([h, 0]);

    var array_test1 = [""];
    var array_test2 = [0];

    for (var i = 0; i < chrom_acum_length.length; i++) {
        var num = i + 1;
        array_test1.push("chr" + num);
        array_test2.push(xScale(chrom_acum_length[i]));
    }
    console.log(x1,x2,y1,y2);
	console.log(data_weight_pvalue);

    //Create scale top
    var xScale_top = d3.scale.ordinal()
        .domain(array_test1)
        .range(array_test2);

    //Define X axis top
    var xAxis_top = d3.svg.axis()
        .scale(xScale_top)
        .orient("top")
        .ticks(0);

    //Define X axis
    var xAxis = d3.svg.axis()
        .scale(xScale)
        .orient("bottom");

    //Define Y axis
    var yAxis = d3.svg.axis()
        .scale(yScale)
        .orient("left")
        .ticks(5);

    // Define tooltip
    var tooltip = d3.select("#containergenome2")
	    .append("div")
	    .style("position", "absolute")
	    .style("z-index", "10")
	    .style("visibility", "hidden")
	    .text("a simple tooltip");

    //Create SVG element
    var svg = d3.select("#containergenome2")
        .append("svg")
        .attr("width", w + margin.left + margin.right)
        .attr("height", h + margin.top + margin.bottom)
        .append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

	// Vertical chromosome lines
    var line_chrom = svg.selectAll("line")
        .data(chrom_acum_length)
        .enter()
        .append("line")
        .attr("class", "linechrom")
        .attr("x1", function(d) {
            return xScale(d);
        })
        .attr("y1", 0)
        .attr("x2", function(d) {
            return xScale(d);
        })
        .attr("y2", h)
        .attr("stroke-width", 1)
        .attr("stroke-dasharray", 5)
        .style("stroke", "black") //stroke-dasharray="5"
        .style("opacity", 0.5);

    //Create circles for the low significance points
    var circle_low = svg.selectAll("circle_low")
        .data(data_low)
        .enter()
        .append("circle")
        .attr("cx", function(d) {
            return xScale(parseFloat(d[0]));
        })
        .attr("cy", function(d) {
            return yScale(parseFloat(d[1]));
        })
        .attr("r", 2)
        .style("fill", function(d) {
            return "#D3D3D3";
        });

	// Create circles for the high significance points
    var circle_high = svg.selectAll("circle_high")
        .data(data_high)
        .enter()
        .append("circle")
        .attr("cx", function(d) {
            return xScale(parseFloat(d[0]));
        })
        .attr("cy", function(d) {
            return yScale(parseFloat(d[1]));
        })
        .attr("r", function(d){
             return d3.max([5,d3.min([10,1.5*parseFloat(d[1])])]);
             })
        .on("mouseover", function(d){return tooltip.style("visibility", "visible").text("gene: " + d[2] + "; " + "log10(Bayes Factor): " + Number((d[1]).toFixed(2)));})
        .on("mousemove", function(){return tooltip.style("top", (event.pageY-10)+"px").style("left",(event.pageX+10)+"px");})
        .on("mouseout", function(){return tooltip.style("visibility", "hidden");})
        .on("click", function(d){
       	window.open("https://biobankengine.stanford.edu/{{namespace}}/awesome?query=" + d[2], '_blank');
    	})
        .style("fill", function(d) {
            return colorScale(parseInt(d[3]));
        });

    //Create X axis
    svg.append("g")
        .attr("class", "manaxis").attr("font-size", "10px")
        .attr("transform", "translate(0," + (h) + ")")
        .call(xAxis);

    svg.append("g").attr("transform", "translate(0," + 0 + ")")
        .attr("class", "xt axis")
        .call(xAxis_top);

    svg.selectAll(".xt text") // select all the text elements for the xaxis
        .attr("transform", function(d) {
            return "translate(" + this.getBBox().height + "," + this.getBBox().height * -0.5 + ")rotate(-45)";
        });

    //Create Y axis
    svg.append("g")
        .attr("class", "manaxis")
        .call(yAxis)
        .append("text")
        .attr("class", "manlabel")
        .attr("transform", "rotate(-90)")
        .attr("x", -10)
        .attr("y", -50)
        .attr("dy", ".71em")
        .style("text-anchor", "end").attr("font-size", "17px")
        .text("log10 Bayes Factor");
});      

})(d3);
