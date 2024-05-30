HTMLWidgets.widget({
  name: 'bipartitenetwork',
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

      var tallestColumnSize = Math.max(x.data.column1NodeIds.length, x.data.column2NodeIds.length)

      function getNodeId(d) {
        return(d);
      }

      function findNodeCX(d, isColumn1) {
        var cx = isColumn1 ? innerWidth / 6 : innerWidth / 6 * 5;
        return(cx)
      }

      function findNodeLabelX(d, isColumn1) {
        var cx = isColumn1 ? findNodeCX(d, isColumn1) - 50 : findNodeCX(d, isColumn1) + 10;
        return(cx)
      }

      function findNodeCY(d, i) {
        return(i * (innerHeight / tallestColumnSize) + 10)
      }

      function findLinkY1(d) {
        const isCurrentSourceNode = (element) => element === d.source;
        var i = x.data.column1NodeIds.findIndex(isCurrentSourceNode);
        var y1 = findNodeCY(d, i);
        return(y1);
      }

      function findLinkY2(d) {
        const isCurrentTargetNode = (element) => element === d.target;
        var i = x.data.column2NodeIds.findIndex(isCurrentTargetNode);
        var y2 = findNodeCY(d, i);
        return(y2);
      }
      
      svg.selectAll('.link')
        .data(x.data.links)
        .enter()
        .append('line')
        .attr('class', 'link')
        .attr('x1', d => findNodeCX(d, true))
        .attr('y1', d => findLinkY1(d))
        .attr('x2', d => findNodeCX(d, false))
        .attr('y2', d => findLinkY2(d))
        .style('stroke', function(d) { return d.value > 0 ? "#0EADA5" : "#AD3C00"; })
        .style('stroke-width', function(d) { return 2 * Math.abs(d.value); });
  
      let sources = 
      svg.selectAll('.node-source')
        .data(x.data.column1NodeIds)
        .enter()
        .append('g')
        .attr('class', 'node');

      sources  
        .append('circle')
        .attr('id', d => getNodeId(d))
        .attr('r', 5)
        .attr('cx', d => findNodeCX(d, true))
        .attr('cy', (d,i) => findNodeCY(d,i))
        .style('fill', 'white')
        .style('stroke', 'black');

      sources  
        .append('text')
        .attr('x', d => findNodeLabelX(d, true))
        .attr('y', (d,i) => findNodeCY(d,i))
        .text(d => getNodeId(d));

      let targets =
      svg.selectAll('.node-target')
        .data(x.data.column2NodeIds)
        .enter()
        .append('g')
        .attr('class', 'node');

      targets
        .append('circle')
        .attr('id', d => getNodeId(d))
        .attr('r', 5)
        .attr('cx', d => findNodeCX(d, false))
        .attr('cy', (d,i) => findNodeCY(d,i))
        .style('fill', 'white')
        .style('stroke', 'black');

      targets
        .append('text')
        .attr('x', d => findNodeLabelX(d, false)) 
        .attr('y', (d,i) => findNodeCY(d,i))
        .text(d => getNodeId(d));
  
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
