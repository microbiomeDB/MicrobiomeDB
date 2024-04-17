HTMLWidgets.widget({
  name: 'unipartitenetwork',
  type: 'output',

  factory: function(el, width, height) {
    var margin = { top: 20, right: 20, bottom: 30, left: 50 },
        innerWidth = width - margin.left - margin.right,
        innerHeight = height - margin.top - margin.bottom;

    var d3 = d3v5;
    
    var svg = d3.select(el).append('svg')
      .attr('width', width)
      .attr('height', height)
      .attr('transform', 'translate(' + margin.left + ',' + margin.top + ')');

    function renderValue(x) {
      svg.selectAll('*').remove();

      // see insipiration here https://observablehq.com/@chakreshiitgn/bipartite-network-of-teams
      var links = x.data.links;
      var nodes = x.data.nodes;

      // Force Diagram
      var simulation = d3.forceSimulation(nodes)
        .force("link", d3.forceLink(links).id(d => d.id))
        .force("charge", d3.forceManyBody().strength(-500))
        .force("center", d3.forceCenter(width / 2, height / 2));

      var link = svg
			  .selectAll(".link")
			  .data(links)
			  .enter()
			  .append("polyline") //Create as polyline for arrows
        .attr('stroke', function(d) { return d.value > 0 ? "#0EADA5" : "#AD3C00"; })
        .attr('stroke-width', function(d) { return 2 * Math.abs(d.value); })
			  .attr("marker-mid", "url(#end)");  // Add Marker

      // Create nodes
      var node = svg
			  .append("g")
			  .selectAll(".node")
			  .data(nodes)
        .enter()
        .append("circle")
        .attr("r", 5)
        .attr('fill', 'white')
        .attr('stroke', 'black')
			  //.attr("r", d=>5+d.degree) // add degree to data later to bring this back
			  .call(d3.drag(simulation)).attr('class', 'node');
  
       var nodeLinkStatus = {};
       links.forEach(d => {
         nodeLinkStatus[`${d.source.index},${d.target.index}`] = 1;
       });

      function isConnected(a, b) {
        return nodeLinkStatus[`${a.index},${b.index}`] || nodeLinkStatus[`${b.index},${a.index}`] || a.index === b.index;
      }

      // Node interactibility
      node.on('mouseover',function (d) {
        node.style('stroke-opacity', function (o) {
          var thisOpacity = 0;
          if(isConnected(d, o)){
            thisOpacity = 1;
          } else{
            thisOpacity = 0.3;
          }
          this.setAttribute('fill-opacity', thisOpacity);

          return thisOpacity;        
        });
    
        link.style('opacity', function(l) {
          if (d === l.source || d === l.target){
            return 1;
          } else{
            return 0.2;
          }
        });
    
   
        var xpos =d.x;
        var ypos = d.y;
        var tgrp = svg.append("g")
          .attr("id", "tooltip")
          .attr("transform", (d, i) => `translate(${xpos+10},${ypos})`);
        tgrp.append("rect")
          .attr("width", "140px")
          .attr("height", "24px")
          .attr("fill", "gray")
        tgrp.append("text")
          .attr("x", 5)
          .attr("y", 14)
          .attr("text-anchor", "left")
          .attr("font-family", "sans-serif")
          .attr("font-size", "11px")
          .attr("font-weight", "bold")
          .attr("fill", "white")
          .text(`${d.id}`);
      });
  
      node.on('mouseout',function (d) {
        node.style('stroke-opacity', function (o) {
          this.setAttribute('fill-opacity', 1);
          return 1;
        });
        link.style('opacity',1);
        link.style('stroke-width', function(d) { return 2 * Math.abs(d.value); });
        d3.select("#tooltip").remove();
      });
   
      // Simulation tick
      simulation.on("tick", () => {
        link.attr("points", function(d) {
          return d.source.x + "," + d.source.y + " " + 
                (d.source.x + d.target.x)/2 + "," + (d.source.y + d.target.y)/2 + " " +
                 d.target.x + "," + d.target.y; });
        
        node
          .attr("cx", d => d.x)
          .attr("cy", d => d.y);
      });

      return svg.node();
    }

    return {
      renderValue: renderValue,

      resize: function(newWidth, newHeight) {
        svg.attr('width', newWidth).attr('height', newHeight);
        innerWidth = newWidth - margin.left - margin.right;
        innerHeight = newHeight - margin.top - margin.bottom;
        // TODO make sure resizing works
        renderValue(svg.datum());
      }
    };
  }
});
