var styles = [ {
  "format_version" : "1.0",
  "generated_by" : "cytoscape-3.7.2",
  "target_cytoscapejs_version" : "~2.1",
  "title" : "BEL Visualization",
  "style" : [ {
    "selector" : "node",
    "css" : {
      "border-width" : 1.0,
      "height" : 30.0,
      "shape" : "rectangle",
      "text-valign" : "center",
      "text-halign" : "center",
      "border-color" : "rgb(0,0,0)",
      "color" : "rgb(0,0,0)",
      "font-size" : 12,
      "font-family" : "Dialog",
      "font-weight" : "normal",
      "width" : 70.0,
      "border-opacity" : 1.0,
      "text-opacity" : 1.0,
      "background-color" : "rgb(255,255,255)",
      "background-opacity" : 1.0,
      "content" : "data(name)"
    }
  }, {
    "selector" : "node[bel_function_type = 'transportActivity']",
    "css" : {
      "width" : 110.0
    }
  }, {
    "selector" : "node[bel_function_type = 'cellSecretion']",
    "css" : {
      "width" : 110.0
    }
  }, {
    "selector" : "node[bel_function_type = 'catalyticActivity']",
    "css" : {
      "width" : 110.0
    }
  }, {
    "selector" : "node[bel_function_type = 'reaction']",
    "css" : {
      "width" : 150.0
    }
  }, {
    "selector" : "node[bel_function_type = 'pathology']",
    "css" : {
      "width" : 125.0
    }
  }, {
    "selector" : "node[bel_function_type = 'gtpBoundActivity']",
    "css" : {
      "width" : 110.0
    }
  }, {
    "selector" : "node[bel_function_type = 'translocation']",
    "css" : {
      "width" : 110.0
    }
  }, {
    "selector" : "node[bel_function_type = 'kinaseActivity']",
    "css" : {
      "width" : 110.0
    }
  }, {
    "selector" : "node[bel_function_type = 'proteinAbundance']",
    "css" : {
      "width" : 75.0
    }
  }, {
    "selector" : "node[bel_function_type = 'peptidaseActivity']",
    "css" : {
      "width" : 110.0
    }
  }, {
    "selector" : "node[bel_function_type = 'abundance']",
    "css" : {
      "width" : 75.0
    }
  }, {
    "selector" : "node[bel_function_type = 'degradation']",
    "css" : {
      "width" : 110.0
    }
  }, {
    "selector" : "node[bel_function_type = 'compositeAbundance']",
    "css" : {
      "width" : 150.0
    }
  }, {
    "selector" : "node[bel_function_type = 'molecularActivity']",
    "css" : {
      "width" : 110.0
    }
  }, {
    "selector" : "node[bel_function_type = 'rnaAbundance']",
    "css" : {
      "width" : 75.0
    }
  }, {
    "selector" : "node[bel_function_type = 'biologicalProcess']",
    "css" : {
      "width" : 125.0
    }
  }, {
    "selector" : "node[bel_function_type = 'microRNAAbundance']",
    "css" : {
      "width" : 75.0
    }
  }, {
    "selector" : "node[bel_function_type = 'geneAbundance']",
    "css" : {
      "width" : 75.0
    }
  }, {
    "selector" : "node[bel_function_type = 'ribosylationActivity']",
    "css" : {
      "width" : 110.0
    }
  }, {
    "selector" : "node[bel_function_type = 'cellSurfaceExpression']",
    "css" : {
      "width" : 110.0
    }
  }, {
    "selector" : "node[bel_function_type = 'chaperoneActivity']",
    "css" : {
      "width" : 110.0
    }
  }, {
    "selector" : "node[bel_function_type = 'transcriptionalActivity']",
    "css" : {
      "width" : 110.0
    }
  }, {
    "selector" : "node[bel_function_type = 'phosphataseActivity']",
    "css" : {
      "width" : 110.0
    }
  }, {
    "selector" : "node[bel_function_type = 'complexAbundance']",
    "css" : {
      "width" : 150.0
    }
  }, {
    "selector" : "node[linked]",
    "css" : {
      "border-opacity" : 1.0
    }
  }, {
    "selector" : "node[!linked]",
    "css" : {
      "border-opacity" : 0.39215686274509803
    }
  }, {
    "selector" : "node[linked]",
    "css" : {
      "text-opacity" : 1.0
    }
  }, {
    "selector" : "node[!linked]",
    "css" : {
      "text-opacity" : 0.39215686274509803
    }
  }, {
    "selector" : "node[bel_function_type = 'transportActivity']",
    "css" : {
      "background-color" : "rgb(100,100,255)"
    }
  }, {
    "selector" : "node[bel_function_type = 'cellSecretion']",
    "css" : {
      "background-color" : "rgb(200,200,255)"
    }
  }, {
    "selector" : "node[bel_function_type = 'catalyticActivity']",
    "css" : {
      "background-color" : "rgb(100,100,255)"
    }
  }, {
    "selector" : "node[bel_function_type = 'reaction']",
    "css" : {
      "background-color" : "rgb(150,150,255)"
    }
  }, {
    "selector" : "node[bel_function_type = 'pathology']",
    "css" : {
      "background-color" : "rgb(255,50,100)"
    }
  }, {
    "selector" : "node[bel_function_type = 'gtpBoundActivity']",
    "css" : {
      "background-color" : "rgb(100,100,255)"
    }
  }, {
    "selector" : "node[bel_function_type = 'translocation']",
    "css" : {
      "background-color" : "rgb(200,200,255)"
    }
  }, {
    "selector" : "node[bel_function_type = 'kinaseActivity']",
    "css" : {
      "background-color" : "rgb(100,100,255)"
    }
  }, {
    "selector" : "node[bel_function_type = 'proteinAbundance']",
    "css" : {
      "background-color" : "rgb(100,255,255)"
    }
  }, {
    "selector" : "node[bel_function_type = 'peptidaseActivity']",
    "css" : {
      "background-color" : "rgb(100,100,255)"
    }
  }, {
    "selector" : "node[bel_function_type = 'abundance']",
    "css" : {
      "background-color" : "rgb(100,200,255)"
    }
  }, {
    "selector" : "node[bel_function_type = 'degradation']",
    "css" : {
      "background-color" : "rgb(255,50,100)"
    }
  }, {
    "selector" : "node[bel_function_type = 'compositeAbundance']",
    "css" : {
      "background-color" : "rgb(200,255,255)"
    }
  }, {
    "selector" : "node[bel_function_type = 'molecularActivity']",
    "css" : {
      "background-color" : "rgb(100,100,255)"
    }
  }, {
    "selector" : "node[bel_function_type = 'rnaAbundance']",
    "css" : {
      "background-color" : "rgb(100,255,100)"
    }
  }, {
    "selector" : "node[bel_function_type = 'biologicalProcess']",
    "css" : {
      "background-color" : "rgb(255,100,255)"
    }
  }, {
    "selector" : "node[bel_function_type = 'microRNAAbundance']",
    "css" : {
      "background-color" : "rgb(150,255,150)"
    }
  }, {
    "selector" : "node[bel_function_type = 'geneAbundance']",
    "css" : {
      "background-color" : "rgb(200,255,200)"
    }
  }, {
    "selector" : "node[bel_function_type = 'ribosylationActivity']",
    "css" : {
      "background-color" : "rgb(100,100,255)"
    }
  }, {
    "selector" : "node[bel_function_type = 'cellSurfaceExpression']",
    "css" : {
      "background-color" : "rgb(200,200,255)"
    }
  }, {
    "selector" : "node[bel_function_type = 'chaperoneActivity']",
    "css" : {
      "background-color" : "rgb(100,100,255)"
    }
  }, {
    "selector" : "node[bel_function_type = 'transcriptionalActivity']",
    "css" : {
      "background-color" : "rgb(100,100,255)"
    }
  }, {
    "selector" : "node[bel_function_type = 'phosphataseActivity']",
    "css" : {
      "background-color" : "rgb(100,100,255)"
    }
  }, {
    "selector" : "node[bel_function_type = 'complexAbundance']",
    "css" : {
      "background-color" : "rgb(100,150,255)"
    }
  }, {
    "selector" : "node[linked]",
    "css" : {
      "background-opacity" : 1.0
    }
  }, {
    "selector" : "node[!linked]",
    "css" : {
      "background-opacity" : 0.39215686274509803
    }
  }, {
    "selector" : "node[bel_function_type = 'transcriptionalActivity']",
    "css" : {
      "width" : 100.0,
      "height" : 100.0
    }
  }, {
    "selector" : "node[bel_function_type = 'transportActivity']",
    "css" : {
      "shape" : "roundrectangle"
    }
  }, {
    "selector" : "node[bel_function_type = 'cellSecretion']",
    "css" : {
      "shape" : "v"
    }
  }, {
    "selector" : "node[bel_function_type = 'catalyticActivity']",
    "css" : {
      "shape" : "hexagon"
    }
  }, {
    "selector" : "node[bel_function_type = 'reaction']",
    "css" : {
      "shape" : "parallelogram"
    }
  }, {
    "selector" : "node[bel_function_type = 'pathology']",
    "css" : {
      "shape" : "parallelogram"
    }
  }, {
    "selector" : "node[bel_function_type = 'gtpBoundActivity']",
    "css" : {
      "shape" : "hexagon"
    }
  }, {
    "selector" : "node[bel_function_type = 'translocation']",
    "css" : {
      "shape" : "v"
    }
  }, {
    "selector" : "node[bel_function_type = 'kinaseActivity']",
    "css" : {
      "shape" : "hexagon"
    }
  }, {
    "selector" : "node[bel_function_type = 'proteinAbundance']",
    "css" : {
      "shape" : "roundrectangle"
    }
  }, {
    "selector" : "node[bel_function_type = 'peptidaseActivity']",
    "css" : {
      "shape" : "hexagon"
    }
  }, {
    "selector" : "node[bel_function_type = 'abundance']",
    "css" : {
      "shape" : "diamond"
    }
  }, {
    "selector" : "node[bel_function_type = 'degradation']",
    "css" : {
      "shape" : "roundrectangle"
    }
  }, {
    "selector" : "node[bel_function_type = 'compositeAbundance']",
    "css" : {
      "shape" : "roundrectangle"
    }
  }, {
    "selector" : "node[bel_function_type = 'molecularActivity']",
    "css" : {
      "shape" : "hexagon"
    }
  }, {
    "selector" : "node[bel_function_type = 'rnaAbundance']",
    "css" : {
      "shape" : "ellipse"
    }
  }, {
    "selector" : "node[bel_function_type = 'biologicalProcess']",
    "css" : {
      "shape" : "parallelogram"
    }
  }, {
    "selector" : "node[bel_function_type = 'microRNAAbundance']",
    "css" : {
      "shape" : "ellipse"
    }
  }, {
    "selector" : "node[bel_function_type = 'geneAbundance']",
    "css" : {
      "shape" : "ellipse"
    }
  }, {
    "selector" : "node[bel_function_type = 'ribosylationActivity']",
    "css" : {
      "shape" : "hexagon"
    }
  }, {
    "selector" : "node[bel_function_type = 'cellSurfaceExpression']",
    "css" : {
      "shape" : "v"
    }
  }, {
    "selector" : "node[bel_function_type = 'chaperoneActivity']",
    "css" : {
      "shape" : "hexagon"
    }
  }, {
    "selector" : "node[bel_function_type = 'transcriptionalActivity']",
    "css" : {
      "shape" : "hexagon"
    }
  }, {
    "selector" : "node[bel_function_type = 'phosphataseActivity']",
    "css" : {
      "shape" : "hexagon"
    }
  }, {
    "selector" : "node[bel_function_type = 'complexAbundance']",
    "css" : {
      "shape" : "roundrectangle"
    }
  }, {
    "selector" : "node[ id = '107' ]",
    "css" : { }
  }, {
    "selector" : "node:selected",
    "css" : {
      "background-color" : "rgb(255,255,0)"
    }
  }, {
    "selector" : "edge",
    "css" : {
      "font-size" : 10,
      "target-arrow-color" : "rgb(0,0,0)",
      "source-arrow-shape" : "none",
      "width" : 1.0,
      "content" : "",
      "opacity" : 1.0,
      "color" : "rgb(0,0,0)",
      "line-color" : "rgb(0,0,0)",
      "line-style" : "solid",
      "target-arrow-shape" : "none",
      "text-opacity" : 1.0,
      "source-arrow-color" : "rgb(0,0,0)",
      "font-family" : "Dialog",
      "font-weight" : "normal"
    }
  }, {
    "selector" : "edge[interaction = 'transcribedTo']",
    "css" : {
      "line-style" : "dashed"
    }
  }, {
    "selector" : "edge[interaction = 'causesNoChange']",
    "css" : {
      "line-style" : "dashed"
    }
  }, {
    "selector" : "edge[interaction = 'translatedTo']",
    "css" : {
      "line-style" : "dashed"
    }
  }, {
    "selector" : "edge[interaction = 'hasProduct']",
    "css" : {
      "line-style" : "dashed"
    }
  }, {
    "selector" : "edge[interaction = 'positiveCorrelation']",
    "css" : {
      "line-style" : "dotted"
    }
  }, {
    "selector" : "edge[interaction = 'reactantIn']",
    "css" : {
      "line-style" : "dashed"
    }
  }, {
    "selector" : "edge[interaction = 'subProcessOf']",
    "css" : {
      "line-style" : "dashed"
    }
  }, {
    "selector" : "edge[interaction = 'hasMember']",
    "css" : {
      "line-style" : "dashed"
    }
  }, {
    "selector" : "edge[interaction = 'orthologous']",
    "css" : {
      "line-style" : "dashed"
    }
  }, {
    "selector" : "edge[interaction = 'association']",
    "css" : {
      "line-style" : "dotted"
    }
  }, {
    "selector" : "edge[interaction = 'actsIn']",
    "css" : {
      "line-style" : "dashed"
    }
  }, {
    "selector" : "edge[interaction = 'includes']",
    "css" : {
      "line-style" : "dashed"
    }
  }, {
    "selector" : "edge[interaction = 'hasVariant']",
    "css" : {
      "line-style" : "dashed"
    }
  }, {
    "selector" : "edge[interaction = 'directlyDecreases']",
    "css" : {
      "line-style" : "solid"
    }
  }, {
    "selector" : "edge[interaction = 'negativeCorrelation']",
    "css" : {
      "line-style" : "dotted"
    }
  }, {
    "selector" : "edge[interaction = 'directlyIncreases']",
    "css" : {
      "line-style" : "solid"
    }
  }, {
    "selector" : "edge[interaction = 'prognosticBiomarkerFor']",
    "css" : {
      "line-style" : "dotted"
    }
  }, {
    "selector" : "edge[interaction = 'hasComponent']",
    "css" : {
      "line-style" : "dashed"
    }
  }, {
    "selector" : "edge[interaction = 'hasModification']",
    "css" : {
      "line-style" : "dashed"
    }
  }, {
    "selector" : "edge[interaction = 'subProcess_OF']",
    "css" : {
      "line-style" : "dashed"
    }
  }, {
    "selector" : "edge[interaction = 'isA']",
    "css" : {
      "line-style" : "dashed"
    }
  }, {
    "selector" : "edge[interaction = 'rateLimitingStepOf']",
    "css" : {
      "line-style" : "dashed"
    }
  }, {
    "selector" : "edge[interaction = 'biomarkerFor']",
    "css" : {
      "line-style" : "dotted"
    }
  }, {
    "selector" : "edge[interaction = 'translocates']",
    "css" : {
      "line-style" : "dashed"
    }
  }, {
    "selector" : "edge[interaction = 'translatedTo']",
    "css" : {
      "target-arrow-shape" : "diamond"
    }
  }, {
    "selector" : "edge[interaction = 'hasProduct']",
    "css" : {
      "target-arrow-shape" : "diamond"
    }
  }, {
    "selector" : "edge[interaction = 'reactantIn']",
    "css" : {
      "target-arrow-shape" : "diamond"
    }
  }, {
    "selector" : "edge[interaction = 'orthologous']",
    "css" : {
      "target-arrow-shape" : "none"
    }
  }, {
    "selector" : "edge[interaction = 'actsIn']",
    "css" : {
      "target-arrow-shape" : "diamond"
    }
  }, {
    "selector" : "edge[interaction = 'association']",
    "css" : {
      "target-arrow-shape" : "none"
    }
  }, {
    "selector" : "edge[interaction = 'hasVariant']",
    "css" : {
      "target-arrow-shape" : "diamond"
    }
  }, {
    "selector" : "edge[interaction = 'directlyIncreases']",
    "css" : {
      "target-arrow-shape" : "triangle"
    }
  }, {
    "selector" : "edge[interaction = 'unknown']",
    "css" : {
      "target-arrow-shape" : "triangle"
    }
  }, {
    "selector" : "edge[interaction = 'increases']",
    "css" : {
      "target-arrow-shape" : "triangle"
    }
  }, {
    "selector" : "edge[interaction = 'prognosticBiomarkerFor']",
    "css" : {
      "target-arrow-shape" : "triangle"
    }
  }, {
    "selector" : "edge[interaction = 'hasComponent']",
    "css" : {
      "target-arrow-shape" : "diamond"
    }
  }, {
    "selector" : "edge[interaction = 'hasModification']",
    "css" : {
      "target-arrow-shape" : "diamond"
    }
  }, {
    "selector" : "edge[interaction = 'subProcess_OF']",
    "css" : {
      "target-arrow-shape" : "diamond"
    }
  }, {
    "selector" : "edge[interaction = 'isA']",
    "css" : {
      "target-arrow-shape" : "diamond"
    }
  }, {
    "selector" : "edge[interaction = 'rateLimitingStepOf']",
    "css" : {
      "target-arrow-shape" : "triangle"
    }
  }, {
    "selector" : "edge[interaction = 'biomarkerFor']",
    "css" : {
      "target-arrow-shape" : "triangle"
    }
  }, {
    "selector" : "edge[interaction = 'translocates']",
    "css" : {
      "target-arrow-shape" : "diamond"
    }
  }, {
    "selector" : "edge[interaction = 'analogous']",
    "css" : {
      "target-arrow-shape" : "triangle"
    }
  }, {
    "selector" : "edge[interaction = 'transcribedTo']",
    "css" : {
      "target-arrow-shape" : "diamond"
    }
  }, {
    "selector" : "edge[interaction = 'positiveCorrelation']",
    "css" : {
      "target-arrow-shape" : "triangle"
    }
  }, {
    "selector" : "edge[interaction = 'subProcessOf']",
    "css" : {
      "target-arrow-shape" : "diamond"
    }
  }, {
    "selector" : "edge[interaction = 'hasMember']",
    "css" : {
      "target-arrow-shape" : "diamond"
    }
  }, {
    "selector" : "edge[interaction = 'includes']",
    "css" : {
      "target-arrow-shape" : "diamond"
    }
  }, {
    "selector" : "edge[interaction = 'directlyDecreases']",
    "css" : {
      "target-arrow-shape" : "circle"
    }
  }, {
    "selector" : "edge[interaction = 'negativeCorrelation']",
    "css" : {
      "target-arrow-shape" : "circle"
    }
  }, {
    "selector" : "edge[interaction = 'decreases']",
    "css" : {
      "target-arrow-shape" : "circle"
    }
  }, {
    "selector" : "edge[interaction = 'positiveCorrelation']",
    "css" : {
      "source-arrow-shape" : "triangle"
    }
  }, {
    "selector" : "edge[interaction = 'hasMember']",
    "css" : {
      "source-arrow-shape" : "none"
    }
  }, {
    "selector" : "edge[interaction = 'orthologous']",
    "css" : {
      "source-arrow-shape" : "none"
    }
  }, {
    "selector" : "edge[interaction = 'association']",
    "css" : {
      "source-arrow-shape" : "none"
    }
  }, {
    "selector" : "edge[interaction = 'negativeCorrelation']",
    "css" : {
      "source-arrow-shape" : "circle"
    }
  }, {
    "selector" : "edge[interaction = 'translatedTo']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'hasProduct']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'reactantIn']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'orthologous']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'actsIn']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'association']",
    "css" : {
      "width" : 1.0
    }
  }, {
    "selector" : "edge[interaction = 'hasVariant']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'directlyIncreases']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'increases']",
    "css" : {
      "width" : 1.0
    }
  }, {
    "selector" : "edge[interaction = 'prognosticBiomarkerFor']",
    "css" : {
      "width" : 1.0
    }
  }, {
    "selector" : "edge[interaction = 'hasComponent']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'hasModification']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'subProcess_OF']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'isA']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'rateLimitingStepOf']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'biomarkerFor']",
    "css" : {
      "width" : 1.0
    }
  }, {
    "selector" : "edge[interaction = 'translocates']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'transcribedTo']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'causesNoChange']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'positiveCorrelation']",
    "css" : {
      "width" : 1.0
    }
  }, {
    "selector" : "edge[interaction = 'subProcessOf']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'hasMember']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'includes']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'directlyDecreases']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'negativeCorrelation']",
    "css" : {
      "width" : 1.0
    }
  }, {
    "selector" : "edge[interaction = 'decreases']",
    "css" : {
      "width" : 1.0
    }
  }, {
    "selector" : "edge[linked]",
    "css" : {
      "opacity" : 1.0
    }
  }, {
    "selector" : "edge[!linked]",
    "css" : {
      "opacity" : 0.19607843137254902
    }
  }, {
    "selector" : "edge[interaction = 'transcribedTo']",
    "css" : {
      "line-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'translatedTo']",
    "css" : {
      "line-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'hasProduct']",
    "css" : {
      "line-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'reactantIn']",
    "css" : {
      "line-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'subProcessOf']",
    "css" : {
      "line-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'hasMember']",
    "css" : {
      "line-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'orthologous']",
    "css" : {
      "line-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'actsIn']",
    "css" : {
      "line-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'includes']",
    "css" : {
      "line-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'hasVariant']",
    "css" : {
      "line-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'hasComponent']",
    "css" : {
      "line-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'hasModification']",
    "css" : {
      "line-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'subProcess_OF']",
    "css" : {
      "line-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'isA']",
    "css" : {
      "line-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'translocates']",
    "css" : {
      "line-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'transcribedTo']",
    "css" : {
      "line-color" : "rgb(153,153,153)",
      "target-arrow-color" : "rgb(153,153,153)",
      "source-arrow-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'translatedTo']",
    "css" : {
      "line-color" : "rgb(153,153,153)",
      "target-arrow-color" : "rgb(153,153,153)",
      "source-arrow-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'hasProduct']",
    "css" : {
      "line-color" : "rgb(153,153,153)",
      "target-arrow-color" : "rgb(153,153,153)",
      "source-arrow-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'reactantIn']",
    "css" : {
      "line-color" : "rgb(153,153,153)",
      "target-arrow-color" : "rgb(153,153,153)",
      "source-arrow-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'subProcessOf']",
    "css" : {
      "line-color" : "rgb(153,153,153)",
      "target-arrow-color" : "rgb(153,153,153)",
      "source-arrow-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'hasMember']",
    "css" : {
      "line-color" : "rgb(153,153,153)",
      "target-arrow-color" : "rgb(153,153,153)",
      "source-arrow-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'orthologous']",
    "css" : {
      "line-color" : "rgb(153,153,153)",
      "target-arrow-color" : "rgb(153,153,153)",
      "source-arrow-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'actsIn']",
    "css" : {
      "line-color" : "rgb(153,153,153)",
      "target-arrow-color" : "rgb(153,153,153)",
      "source-arrow-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'includes']",
    "css" : {
      "line-color" : "rgb(153,153,153)",
      "target-arrow-color" : "rgb(153,153,153)",
      "source-arrow-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'hasVariant']",
    "css" : {
      "line-color" : "rgb(153,153,153)",
      "target-arrow-color" : "rgb(153,153,153)",
      "source-arrow-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'hasComponent']",
    "css" : {
      "line-color" : "rgb(153,153,153)",
      "target-arrow-color" : "rgb(153,153,153)",
      "source-arrow-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'hasModification']",
    "css" : {
      "line-color" : "rgb(153,153,153)",
      "target-arrow-color" : "rgb(153,153,153)",
      "source-arrow-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'subProcess_OF']",
    "css" : {
      "line-color" : "rgb(153,153,153)",
      "target-arrow-color" : "rgb(153,153,153)",
      "source-arrow-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'isA']",
    "css" : {
      "line-color" : "rgb(153,153,153)",
      "target-arrow-color" : "rgb(153,153,153)",
      "source-arrow-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'translocates']",
    "css" : {
      "line-color" : "rgb(153,153,153)",
      "target-arrow-color" : "rgb(153,153,153)",
      "source-arrow-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge:selected",
    "css" : {
      "line-color" : "rgb(255,0,0)"
    }
  } ]
}, {
  "format_version" : "1.0",
  "generated_by" : "cytoscape-3.7.2",
  "target_cytoscapejs_version" : "~2.1",
  "title" : "Gradient1",
  "style" : [ {
    "selector" : "node",
    "css" : {
      "border-width" : 0.0,
      "height" : 30.0,
      "shape" : "ellipse",
      "text-valign" : "bottom",
      "text-halign" : "right",
      "border-color" : "rgb(0,0,0)",
      "color" : "rgb(204,204,204)",
      "font-size" : 8,
      "font-family" : "SansSerif",
      "font-weight" : "normal",
      "width" : 30.0,
      "border-opacity" : 1.0,
      "text-opacity" : 1.0,
      "background-color" : "rgb(0,0,0)",
      "background-opacity" : 1.0,
      "content" : "data(name)"
    }
  }, {
    "selector" : "node[ id = '107' ]",
    "css" : { }
  }, {
    "selector" : "node:selected",
    "css" : {
      "background-color" : "rgb(255,255,0)"
    }
  }, {
    "selector" : "edge",
    "css" : {
      "font-size" : 10,
      "target-arrow-color" : "rgb(0,0,0)",
      "source-arrow-shape" : "none",
      "width" : 1.0,
      "content" : "",
      "opacity" : 1.0,
      "color" : "rgb(0,0,0)",
      "line-color" : "rgb(102,102,102)",
      "line-style" : "solid",
      "target-arrow-shape" : "none",
      "text-opacity" : 1.0,
      "source-arrow-color" : "rgb(0,0,0)",
      "font-family" : "SansSerif",
      "font-weight" : "normal"
    }
  }, {
    "selector" : "edge:selected",
    "css" : {
      "line-color" : "rgb(255,0,0)"
    }
  } ]
}, {
  "format_version" : "1.0",
  "generated_by" : "cytoscape-3.7.2",
  "target_cytoscapejs_version" : "~2.1",
  "title" : "Sample3",
  "style" : [ {
    "selector" : "node",
    "css" : {
      "border-width" : 8.0,
      "height" : 20.0,
      "shape" : "ellipse",
      "text-valign" : "bottom",
      "text-halign" : "right",
      "border-color" : "rgb(255,255,255)",
      "color" : "rgb(206,206,206)",
      "font-size" : 14,
      "font-family" : "SansSerif",
      "font-weight" : "normal",
      "width" : 20.0,
      "border-opacity" : 1.0,
      "text-opacity" : 1.0,
      "background-color" : "rgb(61,154,255)",
      "background-opacity" : 1.0,
      "content" : "data(name)"
    }
  }, {
    "selector" : "node[ id = '107' ]",
    "css" : { }
  }, {
    "selector" : "node:selected",
    "css" : {
      "background-color" : "rgb(255,255,0)"
    }
  }, {
    "selector" : "edge",
    "css" : {
      "font-size" : 10,
      "target-arrow-color" : "rgb(0,0,0)",
      "source-arrow-shape" : "none",
      "width" : 2.0,
      "content" : "",
      "opacity" : 1.0,
      "color" : "rgb(0,0,0)",
      "line-color" : "rgb(255,255,255)",
      "line-style" : "solid",
      "target-arrow-shape" : "none",
      "text-opacity" : 1.0,
      "source-arrow-color" : "rgb(0,0,0)",
      "font-family" : "SansSerif",
      "font-weight" : "normal"
    }
  }, {
    "selector" : "edge:selected",
    "css" : {
      "line-color" : "rgb(255,0,0)"
    }
  } ]
}, {
  "format_version" : "1.0",
  "generated_by" : "cytoscape-3.7.2",
  "target_cytoscapejs_version" : "~2.1",
  "title" : "Universe",
  "style" : [ {
    "selector" : "node",
    "css" : {
      "border-width" : 0.0,
      "height" : 40.0,
      "shape" : "ellipse",
      "text-valign" : "center",
      "text-halign" : "center",
      "border-color" : "rgb(0,0,0)",
      "color" : "rgb(255,255,255)",
      "font-size" : 18,
      "font-family" : "Monospaced",
      "font-weight" : "normal",
      "width" : 40.0,
      "border-opacity" : 1.0,
      "text-opacity" : 1.0,
      "background-color" : "rgb(0,0,0)",
      "background-opacity" : 1.0,
      "content" : "data(name)"
    }
  }, {
    "selector" : "node[ id = '107' ]",
    "css" : { }
  }, {
    "selector" : "node:selected",
    "css" : {
      "background-color" : "rgb(255,255,0)"
    }
  }, {
    "selector" : "edge",
    "css" : {
      "font-size" : 10,
      "target-arrow-color" : "rgb(0,0,0)",
      "source-arrow-shape" : "none",
      "width" : 2.0,
      "content" : "",
      "opacity" : 1.0,
      "color" : "rgb(0,0,0)",
      "line-color" : "rgb(153,153,153)",
      "line-style" : "dashed",
      "target-arrow-shape" : "none",
      "text-opacity" : 1.0,
      "source-arrow-color" : "rgb(0,0,0)",
      "font-family" : "Dialog",
      "font-weight" : "normal"
    }
  }, {
    "selector" : "edge:selected",
    "css" : {
      "line-color" : "rgb(255,0,0)"
    }
  } ]
}, {
  "format_version" : "1.0",
  "generated_by" : "cytoscape-3.7.2",
  "target_cytoscapejs_version" : "~2.1",
  "title" : "Solid",
  "style" : [ {
    "selector" : "node",
    "css" : {
      "border-width" : 0.0,
      "height" : 40.0,
      "shape" : "ellipse",
      "text-valign" : "center",
      "text-halign" : "center",
      "border-color" : "rgb(0,0,0)",
      "color" : "rgb(0,0,0)",
      "font-size" : 14,
      "font-family" : "Dialog",
      "font-weight" : "normal",
      "width" : 40.0,
      "border-opacity" : 1.0,
      "text-opacity" : 1.0,
      "background-color" : "rgb(102,102,102)",
      "background-opacity" : 1.0,
      "content" : "data(name)"
    }
  }, {
    "selector" : "node[ id = '107' ]",
    "css" : { }
  }, {
    "selector" : "node:selected",
    "css" : {
      "background-color" : "rgb(255,255,0)"
    }
  }, {
    "selector" : "edge",
    "css" : {
      "font-size" : 10,
      "target-arrow-color" : "rgb(0,0,0)",
      "source-arrow-shape" : "none",
      "width" : 12.0,
      "opacity" : 1.0,
      "color" : "rgb(0,0,0)",
      "line-color" : "rgb(204,204,204)",
      "line-style" : "solid",
      "target-arrow-shape" : "none",
      "text-opacity" : 1.0,
      "source-arrow-color" : "rgb(0,0,0)",
      "font-family" : "Dialog",
      "font-weight" : "normal",
      "content" : "data(interaction)"
    }
  }, {
    "selector" : "edge:selected",
    "css" : {
      "line-color" : "rgb(255,0,0)"
    }
  } ]
}, {
  "format_version" : "1.0",
  "generated_by" : "cytoscape-3.7.2",
  "target_cytoscapejs_version" : "~2.1",
  "title" : "default black",
  "style" : [ {
    "selector" : "node",
    "css" : {
      "border-width" : 0.0,
      "height" : 15.0,
      "shape" : "ellipse",
      "text-valign" : "bottom",
      "text-halign" : "right",
      "border-color" : "rgb(0,153,0)",
      "color" : "rgb(204,204,204)",
      "font-size" : 12,
      "font-family" : "SansSerif",
      "font-weight" : "normal",
      "width" : 15.0,
      "border-opacity" : 1.0,
      "text-opacity" : 1.0,
      "background-color" : "rgb(255,255,255)",
      "background-opacity" : 1.0,
      "content" : "data(name)"
    }
  }, {
    "selector" : "node[ id = '107' ]",
    "css" : { }
  }, {
    "selector" : "node:selected",
    "css" : {
      "background-color" : "rgb(255,255,0)"
    }
  }, {
    "selector" : "edge",
    "css" : {
      "font-size" : 10,
      "target-arrow-color" : "rgb(0,0,0)",
      "source-arrow-shape" : "none",
      "width" : 2.0,
      "content" : "",
      "opacity" : 1.0,
      "color" : "rgb(0,0,0)",
      "line-color" : "rgb(0,153,0)",
      "line-style" : "solid",
      "target-arrow-shape" : "none",
      "text-opacity" : 1.0,
      "source-arrow-color" : "rgb(0,0,0)",
      "font-family" : "Dialog",
      "font-weight" : "normal"
    }
  }, {
    "selector" : "edge:selected",
    "css" : {
      "line-color" : "rgb(255,0,0)"
    }
  } ]
}, {
  "format_version" : "1.0",
  "generated_by" : "cytoscape-3.7.2",
  "target_cytoscapejs_version" : "~2.1",
  "title" : "default",
  "style" : [ {
    "selector" : "node",
    "css" : {
      "border-width" : 0.0,
      "height" : 45.0,
      "shape" : "roundrectangle",
      "text-valign" : "center",
      "text-halign" : "center",
      "border-color" : "rgb(204,204,204)",
      "color" : "rgb(0,0,0)",
      "font-size" : 12,
      "font-family" : "SansSerif",
      "font-weight" : "normal",
      "width" : 90.0,
      "border-opacity" : 1.0,
      "text-opacity" : 1.0,
      "background-color" : "rgb(137,208,245)",
      "background-opacity" : 1.0,
      "content" : "data(name)"
    }
  }, {
    "selector" : "node[ id = '107' ]",
    "css" : { }
  }, {
    "selector" : "node:selected",
    "css" : {
      "background-color" : "rgb(255,255,0)"
    }
  }, {
    "selector" : "edge",
    "css" : {
      "font-size" : 10,
      "target-arrow-color" : "rgb(0,0,0)",
      "source-arrow-shape" : "none",
      "width" : 2.0,
      "content" : "",
      "opacity" : 1.0,
      "color" : "rgb(0,0,0)",
      "line-color" : "rgb(132,132,132)",
      "line-style" : "solid",
      "target-arrow-shape" : "triangle",
      "text-opacity" : 1.0,
      "source-arrow-color" : "rgb(0,0,0)",
      "font-family" : "Dialog",
      "font-weight" : "normal"
    }
  }, {
    "selector" : "edge[is_skeleton_edge]",
    "css" : {
      "line-style" : "solid"
    }
  }, {
    "selector" : "edge[!is_skeleton_edge]",
    "css" : {
      "line-style" : "dotted"
    }
  }, {
    "selector" : "edge[direction_mrpresso = 'Down']",
    "css" : {
      "line-color" : "rgb(255,0,0)"
    }
  }, {
    "selector" : "edge[direction_mrpresso = 'Up']",
    "css" : {
      "line-color" : "rgb(0,0,255)"
    }
  }, {
    "selector" : "edge:selected",
    "css" : {
      "line-color" : "rgb(255,0,0)"
    }
  } ]
}, {
  "format_version" : "1.0",
  "generated_by" : "cytoscape-3.7.2",
  "target_cytoscapejs_version" : "~2.1",
  "title" : "BioPAX_SIF",
  "style" : [ {
    "selector" : "node",
    "css" : {
      "border-width" : 2.0,
      "height" : 40.0,
      "shape" : "ellipse",
      "text-valign" : "center",
      "text-halign" : "center",
      "border-color" : "rgb(0,0,0)",
      "color" : "rgb(0,0,0)",
      "font-size" : 12,
      "font-family" : "SansSerif",
      "font-weight" : "normal",
      "width" : 60.0,
      "border-opacity" : 1.0,
      "text-opacity" : 1.0,
      "background-color" : "rgb(255,153,153)",
      "background-opacity" : 0.49019607843137253,
      "content" : "data(name)"
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'Complex']",
    "css" : {
      "shape" : "hexagon"
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'Complex']",
    "css" : {
      "background-color" : "rgb(153,204,255)"
    }
  }, {
    "selector" : "node[ id = '107' ]",
    "css" : { }
  }, {
    "selector" : "node:selected",
    "css" : {
      "background-color" : "rgb(255,255,0)"
    }
  }, {
    "selector" : "edge",
    "css" : {
      "font-size" : 10,
      "target-arrow-color" : "rgb(0,0,0)",
      "source-arrow-shape" : "none",
      "width" : 4.0,
      "content" : "",
      "opacity" : 1.0,
      "color" : "rgb(0,0,0)",
      "line-color" : "rgb(0,0,0)",
      "line-style" : "solid",
      "target-arrow-shape" : "none",
      "text-opacity" : 1.0,
      "source-arrow-color" : "rgb(0,0,0)",
      "font-family" : "SansSerif",
      "font-weight" : "normal"
    }
  }, {
    "selector" : "edge[interaction = 'controls-expression-of']",
    "css" : {
      "target-arrow-shape" : "triangle"
    }
  }, {
    "selector" : "edge[interaction = 'chemical-affects']",
    "css" : {
      "target-arrow-shape" : "triangle"
    }
  }, {
    "selector" : "edge[interaction = 'controls-state-change-of']",
    "css" : {
      "target-arrow-shape" : "triangle"
    }
  }, {
    "selector" : "edge[interaction = 'controls-phosphorylation-of']",
    "css" : {
      "target-arrow-shape" : "triangle"
    }
  }, {
    "selector" : "edge[interaction = 'controls-transport-of-chemical']",
    "css" : {
      "target-arrow-shape" : "triangle"
    }
  }, {
    "selector" : "edge[interaction = 'used-to-produce']",
    "css" : {
      "target-arrow-shape" : "triangle"
    }
  }, {
    "selector" : "edge[interaction = 'controls-transport-of']",
    "css" : {
      "target-arrow-shape" : "triangle"
    }
  }, {
    "selector" : "edge[interaction = 'consumption-controled-by']",
    "css" : {
      "target-arrow-shape" : "triangle"
    }
  }, {
    "selector" : "edge[interaction = 'controls-production-of']",
    "css" : {
      "target-arrow-shape" : "triangle"
    }
  }, {
    "selector" : "edge[interaction = 'catalysis-precedes']",
    "css" : {
      "target-arrow-shape" : "triangle"
    }
  }, {
    "selector" : "edge[interaction = 'interacts-with']",
    "css" : {
      "line-color" : "rgb(0,85,0)",
      "target-arrow-color" : "rgb(0,85,0)",
      "source-arrow-color" : "rgb(0,85,0)"
    }
  }, {
    "selector" : "edge[interaction = 'chemical-affects']",
    "css" : {
      "line-color" : "rgb(240,144,0)",
      "target-arrow-color" : "rgb(240,144,0)",
      "source-arrow-color" : "rgb(240,144,0)"
    }
  }, {
    "selector" : "edge[interaction = 'controls-state-change-of']",
    "css" : {
      "line-color" : "rgb(0,0,192)",
      "target-arrow-color" : "rgb(0,0,192)",
      "source-arrow-color" : "rgb(0,0,192)"
    }
  }, {
    "selector" : "edge[interaction = 'controls-transport-of']",
    "css" : {
      "line-color" : "rgb(112,0,0)",
      "target-arrow-color" : "rgb(112,0,0)",
      "source-arrow-color" : "rgb(112,0,0)"
    }
  }, {
    "selector" : "edge[interaction = 'consumption-controled-by']",
    "css" : {
      "line-color" : "rgb(255,51,0)",
      "target-arrow-color" : "rgb(255,51,0)",
      "source-arrow-color" : "rgb(255,51,0)"
    }
  }, {
    "selector" : "edge[interaction = 'reacts-with']",
    "css" : {
      "line-color" : "rgb(0,255,0)",
      "target-arrow-color" : "rgb(0,255,0)",
      "source-arrow-color" : "rgb(0,255,0)"
    }
  }, {
    "selector" : "edge[interaction = 'neighbor-of']",
    "css" : {
      "line-color" : "rgb(0,170,0)",
      "target-arrow-color" : "rgb(0,170,0)",
      "source-arrow-color" : "rgb(0,170,0)"
    }
  }, {
    "selector" : "edge[interaction = 'controls-expression-of']",
    "css" : {
      "line-color" : "rgb(0,160,160)",
      "target-arrow-color" : "rgb(0,160,160)",
      "source-arrow-color" : "rgb(0,160,160)"
    }
  }, {
    "selector" : "edge[interaction = 'controls-phosphorylation-of']",
    "css" : {
      "line-color" : "rgb(0,0,255)",
      "target-arrow-color" : "rgb(0,0,255)",
      "source-arrow-color" : "rgb(0,0,255)"
    }
  }, {
    "selector" : "edge[interaction = 'controls-transport-of-chemical']",
    "css" : {
      "line-color" : "rgb(160,0,0)",
      "target-arrow-color" : "rgb(160,0,0)",
      "source-arrow-color" : "rgb(160,0,0)"
    }
  }, {
    "selector" : "edge[interaction = 'used-to-produce']",
    "css" : {
      "line-color" : "rgb(247,85,0)",
      "target-arrow-color" : "rgb(247,85,0)",
      "source-arrow-color" : "rgb(247,85,0)"
    }
  }, {
    "selector" : "edge[interaction = 'in-complex-with']",
    "css" : {
      "line-color" : "rgb(240,0,160)",
      "target-arrow-color" : "rgb(240,0,160)",
      "source-arrow-color" : "rgb(240,0,160)"
    }
  }, {
    "selector" : "edge[interaction = 'controls-production-of']",
    "css" : {
      "line-color" : "rgb(0,204,240)",
      "target-arrow-color" : "rgb(0,204,240)",
      "source-arrow-color" : "rgb(0,204,240)"
    }
  }, {
    "selector" : "edge[interaction = 'catalysis-precedes']",
    "css" : {
      "line-color" : "rgb(112,0,160)",
      "target-arrow-color" : "rgb(112,0,160)",
      "source-arrow-color" : "rgb(112,0,160)"
    }
  }, {
    "selector" : "edge[interaction = 'interacts-with']",
    "css" : {
      "line-color" : "rgb(0,85,0)"
    }
  }, {
    "selector" : "edge[interaction = 'chemical-affects']",
    "css" : {
      "line-color" : "rgb(240,144,0)"
    }
  }, {
    "selector" : "edge[interaction = 'controls-state-change-of']",
    "css" : {
      "line-color" : "rgb(0,0,192)"
    }
  }, {
    "selector" : "edge[interaction = 'controls-transport-of']",
    "css" : {
      "line-color" : "rgb(112,0,0)"
    }
  }, {
    "selector" : "edge[interaction = 'consumption-controled-by']",
    "css" : {
      "line-color" : "rgb(255,51,0)"
    }
  }, {
    "selector" : "edge[interaction = 'reacts-with']",
    "css" : {
      "line-color" : "rgb(0,255,0)"
    }
  }, {
    "selector" : "edge[interaction = 'neighbor-of']",
    "css" : {
      "line-color" : "rgb(0,170,0)"
    }
  }, {
    "selector" : "edge[interaction = 'controls-expression-of']",
    "css" : {
      "line-color" : "rgb(0,160,160)"
    }
  }, {
    "selector" : "edge[interaction = 'controls-phosphorylation-of']",
    "css" : {
      "line-color" : "rgb(0,0,255)"
    }
  }, {
    "selector" : "edge[interaction = 'controls-transport-of-chemical']",
    "css" : {
      "line-color" : "rgb(160,0,0)"
    }
  }, {
    "selector" : "edge[interaction = 'used-to-produce']",
    "css" : {
      "line-color" : "rgb(247,85,0)"
    }
  }, {
    "selector" : "edge[interaction = 'in-complex-with']",
    "css" : {
      "line-color" : "rgb(240,0,160)"
    }
  }, {
    "selector" : "edge[interaction = 'controls-production-of']",
    "css" : {
      "line-color" : "rgb(0,204,240)"
    }
  }, {
    "selector" : "edge[interaction = 'catalysis-precedes']",
    "css" : {
      "line-color" : "rgb(112,0,160)"
    }
  }, {
    "selector" : "edge:selected",
    "css" : {
      "line-color" : "rgb(255,0,0)"
    }
  } ]
}, {
  "format_version" : "1.0",
  "generated_by" : "cytoscape-3.7.2",
  "target_cytoscapejs_version" : "~2.1",
  "title" : "Minimal",
  "style" : [ {
    "selector" : "node",
    "css" : {
      "border-width" : 0.0,
      "height" : 42.0,
      "shape" : "rectangle",
      "text-valign" : "center",
      "text-halign" : "center",
      "border-color" : "rgb(0,0,0)",
      "color" : "rgb(51,51,51)",
      "font-size" : 9,
      "font-family" : "SansSerif",
      "font-weight" : "normal",
      "width" : 42.0,
      "border-opacity" : 1.0,
      "text-opacity" : 1.0,
      "background-color" : "rgb(255,255,255)",
      "background-opacity" : 1.0,
      "content" : "data(name)"
    }
  }, {
    "selector" : "node[ id = '107' ]",
    "css" : { }
  }, {
    "selector" : "node:selected",
    "css" : {
      "background-color" : "rgb(255,255,0)"
    }
  }, {
    "selector" : "edge",
    "css" : {
      "font-size" : 10,
      "target-arrow-color" : "rgb(0,0,0)",
      "source-arrow-shape" : "none",
      "width" : 2.0,
      "content" : "",
      "opacity" : 1.0,
      "color" : "rgb(0,0,0)",
      "line-color" : "rgb(76,76,76)",
      "line-style" : "solid",
      "target-arrow-shape" : "none",
      "text-opacity" : 1.0,
      "source-arrow-color" : "rgb(0,0,0)",
      "font-family" : "SansSerif",
      "font-weight" : "normal"
    }
  }, {
    "selector" : "edge:selected",
    "css" : {
      "line-color" : "rgb(255,0,0)"
    }
  } ]
}, {
  "format_version" : "1.0",
  "generated_by" : "cytoscape-3.7.2",
  "target_cytoscapejs_version" : "~2.1",
  "title" : "Curved",
  "style" : [ {
    "selector" : "node",
    "css" : {
      "border-width" : 7.0,
      "height" : 18.0,
      "shape" : "ellipse",
      "text-valign" : "bottom",
      "text-halign" : "right",
      "border-color" : "rgb(255,255,255)",
      "color" : "rgb(102,102,102)",
      "font-size" : 14,
      "font-family" : "SansSerif",
      "font-weight" : "normal",
      "width" : 18.0,
      "border-opacity" : 1.0,
      "text-opacity" : 1.0,
      "background-color" : "rgb(254,196,79)",
      "background-opacity" : 1.0,
      "content" : "data(name)"
    }
  }, {
    "selector" : "node[ id = '107' ]",
    "css" : { }
  }, {
    "selector" : "node:selected",
    "css" : {
      "background-color" : "rgb(255,255,0)"
    }
  }, {
    "selector" : "edge",
    "css" : {
      "font-size" : 10,
      "target-arrow-color" : "rgb(255,255,255)",
      "source-arrow-shape" : "none",
      "width" : 3.0,
      "content" : "",
      "opacity" : 1.0,
      "color" : "rgb(0,0,0)",
      "line-color" : "rgb(255,255,255)",
      "line-style" : "solid",
      "target-arrow-shape" : "triangle",
      "text-opacity" : 1.0,
      "source-arrow-color" : "rgb(255,255,255)",
      "font-family" : "SansSerif",
      "font-weight" : "normal"
    }
  }, {
    "selector" : "edge:selected",
    "css" : {
      "line-color" : "rgb(255,0,0)"
    }
  } ]
}, {
  "format_version" : "1.0",
  "generated_by" : "cytoscape-3.7.2",
  "target_cytoscapejs_version" : "~2.1",
  "title" : "BioPAX",
  "style" : [ {
    "selector" : "node",
    "css" : {
      "border-width" : 2.0,
      "height" : 20.0,
      "shape" : "ellipse",
      "text-valign" : "center",
      "text-halign" : "center",
      "border-color" : "rgb(0,102,102)",
      "color" : "rgb(0,0,0)",
      "font-size" : 12,
      "font-family" : "SansSerif",
      "font-weight" : "normal",
      "width" : 20.0,
      "border-opacity" : 1.0,
      "text-opacity" : 1.0,
      "background-color" : "rgb(255,255,255)",
      "background-opacity" : 1.0,
      "content" : "data(name)"
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'BiochemicalReaction']",
    "css" : {
      "width" : 13.4
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'GeneticInteraction']",
    "css" : {
      "width" : 13.4
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'Interaction']",
    "css" : {
      "width" : 13.4
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'TransportWithBiochemicalReaction']",
    "css" : {
      "width" : 13.4
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'Conversion']",
    "css" : {
      "width" : 13.4
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'ComplexAssembly']",
    "css" : {
      "width" : 13.4
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'Complex']",
    "css" : {
      "width" : 13.4
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'Degradation']",
    "css" : {
      "width" : 13.4
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'Control']",
    "css" : {
      "width" : 13.4
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'TemplateReactionRegulation']",
    "css" : {
      "width" : 13.4
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'Modulation']",
    "css" : {
      "width" : 13.4
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'MolecularInteraction']",
    "css" : {
      "width" : 13.4
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'TemplateReaction']",
    "css" : {
      "width" : 13.4
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'Catalysis']",
    "css" : {
      "width" : 13.4
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'Transport']",
    "css" : {
      "width" : 13.4
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'SimplePhysicalEntity']",
    "css" : {
      "shape" : "ellipse"
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'Rna']",
    "css" : {
      "shape" : "ellipse"
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'BiochemicalReaction']",
    "css" : {
      "shape" : "rectangle"
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'GeneticInteraction']",
    "css" : {
      "shape" : "rectangle"
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'Interaction']",
    "css" : {
      "shape" : "rectangle"
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'TransportWithBiochemicalReaction']",
    "css" : {
      "shape" : "rectangle"
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'Conversion']",
    "css" : {
      "shape" : "rectangle"
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'ComplexAssembly']",
    "css" : {
      "shape" : "rectangle"
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'Protein']",
    "css" : {
      "shape" : "ellipse"
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'Complex']",
    "css" : {
      "shape" : "diamond"
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'RnaRegion']",
    "css" : {
      "shape" : "ellipse"
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'Degradation']",
    "css" : {
      "shape" : "rectangle"
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'Control']",
    "css" : {
      "shape" : "triangle"
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'TemplateReactionRegulation']",
    "css" : {
      "shape" : "triangle"
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'DnaRegion']",
    "css" : {
      "shape" : "ellipse"
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'PhysicalEntity']",
    "css" : {
      "shape" : "ellipse"
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'SmallMolecule']",
    "css" : {
      "shape" : "ellipse"
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'Dna']",
    "css" : {
      "shape" : "ellipse"
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'Modulation']",
    "css" : {
      "shape" : "triangle"
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'MolecularInteraction']",
    "css" : {
      "shape" : "rectangle"
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'TemplateReaction']",
    "css" : {
      "shape" : "rectangle"
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'Catalysis']",
    "css" : {
      "shape" : "triangle"
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'Transport']",
    "css" : {
      "shape" : "rectangle"
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'Protein-phosphorylated']",
    "css" : {
      "shape" : "ellipse"
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'Complex']",
    "css" : {
      "border-color" : "rgb(0,102,102)"
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'Complex']",
    "css" : {
      "background-color" : "rgb(255,255,255)"
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'BiochemicalReaction']",
    "css" : {
      "height" : 13.4
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'GeneticInteraction']",
    "css" : {
      "height" : 13.4
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'Interaction']",
    "css" : {
      "height" : 13.4
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'TransportWithBiochemicalReaction']",
    "css" : {
      "height" : 13.4
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'Conversion']",
    "css" : {
      "height" : 13.4
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'ComplexAssembly']",
    "css" : {
      "height" : 13.4
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'Complex']",
    "css" : {
      "height" : 13.4
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'Degradation']",
    "css" : {
      "height" : 13.4
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'Control']",
    "css" : {
      "height" : 13.4
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'TemplateReactionRegulation']",
    "css" : {
      "height" : 13.4
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'Modulation']",
    "css" : {
      "height" : 13.4
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'MolecularInteraction']",
    "css" : {
      "height" : 13.4
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'TemplateReaction']",
    "css" : {
      "height" : 13.4
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'Catalysis']",
    "css" : {
      "height" : 13.4
    }
  }, {
    "selector" : "node[BIOPAX_TYPE = 'Transport']",
    "css" : {
      "height" : 13.4
    }
  }, {
    "selector" : "node[ id = '107' ]",
    "css" : { }
  }, {
    "selector" : "node:selected",
    "css" : {
      "background-color" : "rgb(255,255,0)"
    }
  }, {
    "selector" : "edge",
    "css" : {
      "font-size" : 10,
      "target-arrow-color" : "rgb(64,64,64)",
      "source-arrow-shape" : "none",
      "width" : 1.0,
      "content" : "",
      "opacity" : 1.0,
      "color" : "rgb(0,0,0)",
      "line-color" : "rgb(64,64,64)",
      "line-style" : "solid",
      "target-arrow-shape" : "none",
      "text-opacity" : 1.0,
      "source-arrow-color" : "rgb(64,64,64)",
      "font-family" : "SansSerif",
      "font-weight" : "normal"
    }
  }, {
    "selector" : "edge[interaction = 'INHIBITION_NONCOMPETITIVE']",
    "css" : {
      "target-arrow-shape" : "tee"
    }
  }, {
    "selector" : "edge[interaction = 'INHIBITION_OTHER']",
    "css" : {
      "target-arrow-shape" : "tee"
    }
  }, {
    "selector" : "edge[interaction = 'ACTIVATION']",
    "css" : {
      "target-arrow-shape" : "triangle"
    }
  }, {
    "selector" : "edge[interaction = 'INHIBITION_UNCOMPETITIVE']",
    "css" : {
      "target-arrow-shape" : "tee"
    }
  }, {
    "selector" : "edge[interaction = 'cofactor']",
    "css" : {
      "target-arrow-shape" : "triangle"
    }
  }, {
    "selector" : "edge[interaction = 'ACTIVATION_ALLOSTERIC']",
    "css" : {
      "target-arrow-shape" : "triangle"
    }
  }, {
    "selector" : "edge[interaction = 'right']",
    "css" : {
      "target-arrow-shape" : "triangle"
    }
  }, {
    "selector" : "edge[interaction = 'INHIBITION_ALLOSTERIC']",
    "css" : {
      "target-arrow-shape" : "tee"
    }
  }, {
    "selector" : "edge[interaction = 'controlled']",
    "css" : {
      "target-arrow-shape" : "triangle"
    }
  }, {
    "selector" : "edge[interaction = 'contains']",
    "css" : {
      "target-arrow-shape" : "circle"
    }
  }, {
    "selector" : "edge[interaction = 'INHIBITION']",
    "css" : {
      "target-arrow-shape" : "tee"
    }
  }, {
    "selector" : "edge[interaction = 'INHIBITION_UNKMECH']",
    "css" : {
      "target-arrow-shape" : "tee"
    }
  }, {
    "selector" : "edge[interaction = 'INHIBITION_IRREVERSIBLE']",
    "css" : {
      "target-arrow-shape" : "tee"
    }
  }, {
    "selector" : "edge[interaction = 'INHIBITION_COMPETITIVE']",
    "css" : {
      "target-arrow-shape" : "tee"
    }
  }, {
    "selector" : "edge[interaction = 'ACTIVATION_UNKMECH']",
    "css" : {
      "target-arrow-shape" : "triangle"
    }
  }, {
    "selector" : "edge[interaction = 'ACTIVATION_NONALLOSTERIC']",
    "css" : {
      "target-arrow-shape" : "triangle"
    }
  }, {
    "selector" : "edge:selected",
    "css" : {
      "line-color" : "rgb(255,0,0)"
    }
  } ]
}, {
  "format_version" : "1.0",
  "generated_by" : "cytoscape-3.7.2",
  "target_cytoscapejs_version" : "~2.1",
  "title" : "Ripple",
  "style" : [ {
    "selector" : "node",
    "css" : {
      "border-width" : 20.0,
      "height" : 50.0,
      "shape" : "ellipse",
      "text-valign" : "center",
      "text-halign" : "center",
      "border-color" : "rgb(51,153,255)",
      "color" : "rgb(19,58,96)",
      "font-size" : 8,
      "font-family" : "SansSerif",
      "font-weight" : "normal",
      "width" : 50.0,
      "border-opacity" : 1.0,
      "text-opacity" : 1.0,
      "background-color" : "rgb(255,255,255)",
      "background-opacity" : 1.0,
      "content" : "data(name)"
    }
  }, {
    "selector" : "node[ id = '107' ]",
    "css" : { }
  }, {
    "selector" : "node:selected",
    "css" : {
      "background-color" : "rgb(255,255,204)"
    }
  }, {
    "selector" : "edge",
    "css" : {
      "font-size" : 10,
      "target-arrow-color" : "rgb(0,0,0)",
      "source-arrow-shape" : "none",
      "width" : 3.0,
      "content" : "",
      "opacity" : 1.0,
      "color" : "rgb(0,0,0)",
      "line-color" : "rgb(51,153,255)",
      "line-style" : "solid",
      "target-arrow-shape" : "none",
      "text-opacity" : 1.0,
      "source-arrow-color" : "rgb(0,0,0)",
      "font-family" : "SansSerif",
      "font-weight" : "normal"
    }
  }, {
    "selector" : "edge:selected",
    "css" : {
      "line-color" : "rgb(255,0,0)"
    }
  } ]
}, {
  "format_version" : "1.0",
  "generated_by" : "cytoscape-3.7.2",
  "target_cytoscapejs_version" : "~2.1",
  "title" : "size_rank",
  "style" : [ {
    "selector" : "node",
    "css" : {
      "border-width" : 0.0,
      "height" : 12.0,
      "shape" : "rectangle",
      "text-valign" : "center",
      "text-halign" : "center",
      "border-color" : "rgb(0,0,0)",
      "color" : "rgb(51,51,51)",
      "font-size" : 9,
      "font-family" : "SansSerif",
      "font-weight" : "normal",
      "width" : 12.0,
      "border-opacity" : 1.0,
      "text-opacity" : 1.0,
      "background-color" : "rgb(204,204,255)",
      "background-opacity" : 1.0,
      "content" : "data(name)"
    }
  }, {
    "selector" : "node[ id = '107' ]",
    "css" : { }
  }, {
    "selector" : "node:selected",
    "css" : {
      "background-color" : "rgb(255,255,0)"
    }
  }, {
    "selector" : "edge",
    "css" : {
      "font-size" : 10,
      "target-arrow-color" : "rgb(0,0,0)",
      "source-arrow-shape" : "none",
      "width" : 2.0,
      "content" : "",
      "opacity" : 1.0,
      "color" : "rgb(0,0,0)",
      "line-color" : "rgb(76,76,76)",
      "line-style" : "solid",
      "target-arrow-shape" : "none",
      "text-opacity" : 1.0,
      "source-arrow-color" : "rgb(0,0,0)",
      "font-family" : "SansSerif",
      "font-weight" : "normal"
    }
  }, {
    "selector" : "edge:selected",
    "css" : {
      "line-color" : "rgb(255,0,0)"
    }
  } ]
}, {
  "format_version" : "1.0",
  "generated_by" : "cytoscape-3.7.2",
  "target_cytoscapejs_version" : "~2.1",
  "title" : "Big Labels",
  "style" : [ {
    "selector" : "node",
    "css" : {
      "border-width" : 0.0,
      "height" : 5.0,
      "shape" : "ellipse",
      "text-valign" : "center",
      "text-halign" : "center",
      "border-color" : "rgb(0,0,0)",
      "color" : "rgb(51,51,51)",
      "font-size" : 24,
      "font-family" : "SansSerif",
      "font-weight" : "normal",
      "width" : 5.0,
      "border-opacity" : 1.0,
      "text-opacity" : 1.0,
      "background-color" : "rgb(255,255,255)",
      "background-opacity" : 1.0,
      "content" : "data(name)"
    }
  }, {
    "selector" : "node[ id = '107' ]",
    "css" : { }
  }, {
    "selector" : "node:selected",
    "css" : {
      "background-color" : "rgb(255,0,102)"
    }
  }, {
    "selector" : "edge",
    "css" : {
      "font-size" : 10,
      "target-arrow-color" : "rgb(0,0,0)",
      "source-arrow-shape" : "none",
      "width" : 1.0,
      "content" : "",
      "opacity" : 1.0,
      "color" : "rgb(0,0,0)",
      "line-color" : "rgb(183,183,183)",
      "line-style" : "solid",
      "target-arrow-shape" : "none",
      "text-opacity" : 1.0,
      "source-arrow-color" : "rgb(0,0,0)",
      "font-family" : "SansSerif",
      "font-weight" : "normal"
    }
  }, {
    "selector" : "edge:selected",
    "css" : {
      "line-color" : "rgb(255,0,0)"
    }
  } ]
}, {
  "format_version" : "1.0",
  "generated_by" : "cytoscape-3.7.2",
  "target_cytoscapejs_version" : "~2.1",
  "title" : "Nested Network Style",
  "style" : [ {
    "selector" : "node",
    "css" : {
      "border-width" : 2.0,
      "height" : 40.0,
      "shape" : "ellipse",
      "text-valign" : "center",
      "text-halign" : "center",
      "border-color" : "rgb(0,0,0)",
      "color" : "rgb(0,0,0)",
      "font-size" : 12,
      "font-family" : "SansSerif",
      "font-weight" : "normal",
      "width" : 60.0,
      "border-opacity" : 1.0,
      "text-opacity" : 1.0,
      "background-color" : "rgb(255,255,255)",
      "background-opacity" : 1.0,
      "content" : "data(shared_name)"
    }
  }, {
    "selector" : "node[has_nested_network]",
    "css" : {
      "shape" : "rectangle"
    }
  }, {
    "selector" : "node[has_nested_network]",
    "css" : {
      "border-color" : "rgb(0,102,204)"
    }
  }, {
    "selector" : "node[has_nested_network]",
    "css" : {
      "text-valign" : "bottom"
    }
  }, {
    "selector" : "node[has_nested_network]",
    "css" : {
      "background-color" : "rgb(255,255,255)"
    }
  }, {
    "selector" : "node[has_nested_network]",
    "css" : {
      "color" : "rgb(0,102,204)"
    }
  }, {
    "selector" : "node[ id = '107' ]",
    "css" : { }
  }, {
    "selector" : "node:selected",
    "css" : {
      "background-color" : "rgb(255,255,0)"
    }
  }, {
    "selector" : "edge",
    "css" : {
      "font-size" : 10,
      "target-arrow-color" : "rgb(0,0,0)",
      "source-arrow-shape" : "none",
      "width" : 1.0,
      "content" : "",
      "opacity" : 1.0,
      "color" : "rgb(0,0,0)",
      "line-color" : "rgb(64,64,64)",
      "line-style" : "solid",
      "target-arrow-shape" : "none",
      "text-opacity" : 1.0,
      "source-arrow-color" : "rgb(0,0,0)",
      "font-family" : "SansSerif",
      "font-weight" : "normal"
    }
  }, {
    "selector" : "edge:selected",
    "css" : {
      "line-color" : "rgb(255,0,0)"
    }
  } ]
}, {
  "format_version" : "1.0",
  "generated_by" : "cytoscape-3.7.2",
  "target_cytoscapejs_version" : "~2.1",
  "title" : "BEL Visualization Minimal",
  "style" : [ {
    "selector" : "node",
    "css" : {
      "border-width" : 1.0,
      "height" : 30.0,
      "shape" : "rectangle",
      "text-valign" : "center",
      "text-halign" : "center",
      "border-color" : "rgb(0,0,0)",
      "color" : "rgb(0,0,0)",
      "font-size" : 12,
      "font-family" : "Dialog",
      "font-weight" : "normal",
      "width" : 70.0,
      "border-opacity" : 1.0,
      "text-opacity" : 1.0,
      "background-color" : "rgb(255,255,255)",
      "background-opacity" : 1.0,
      "content" : "data(name)"
    }
  }, {
    "selector" : "node[bel_function_type = 'transportActivity']",
    "css" : {
      "width" : 110.0
    }
  }, {
    "selector" : "node[bel_function_type = 'cellSecretion']",
    "css" : {
      "width" : 110.0
    }
  }, {
    "selector" : "node[bel_function_type = 'catalyticActivity']",
    "css" : {
      "width" : 110.0
    }
  }, {
    "selector" : "node[bel_function_type = 'reaction']",
    "css" : {
      "width" : 150.0
    }
  }, {
    "selector" : "node[bel_function_type = 'pathology']",
    "css" : {
      "width" : 125.0
    }
  }, {
    "selector" : "node[bel_function_type = 'gtpBoundActivity']",
    "css" : {
      "width" : 110.0
    }
  }, {
    "selector" : "node[bel_function_type = 'translocation']",
    "css" : {
      "width" : 110.0
    }
  }, {
    "selector" : "node[bel_function_type = 'kinaseActivity']",
    "css" : {
      "width" : 110.0
    }
  }, {
    "selector" : "node[bel_function_type = 'proteinAbundance']",
    "css" : {
      "width" : 75.0
    }
  }, {
    "selector" : "node[bel_function_type = 'peptidaseActivity']",
    "css" : {
      "width" : 110.0
    }
  }, {
    "selector" : "node[bel_function_type = 'abundance']",
    "css" : {
      "width" : 75.0
    }
  }, {
    "selector" : "node[bel_function_type = 'degradation']",
    "css" : {
      "width" : 110.0
    }
  }, {
    "selector" : "node[bel_function_type = 'compositeAbundance']",
    "css" : {
      "width" : 150.0
    }
  }, {
    "selector" : "node[bel_function_type = 'molecularActivity']",
    "css" : {
      "width" : 110.0
    }
  }, {
    "selector" : "node[bel_function_type = 'rnaAbundance']",
    "css" : {
      "width" : 75.0
    }
  }, {
    "selector" : "node[bel_function_type = 'biologicalProcess']",
    "css" : {
      "width" : 125.0
    }
  }, {
    "selector" : "node[bel_function_type = 'microRNAAbundance']",
    "css" : {
      "width" : 75.0
    }
  }, {
    "selector" : "node[bel_function_type = 'geneAbundance']",
    "css" : {
      "width" : 75.0
    }
  }, {
    "selector" : "node[bel_function_type = 'ribosylationActivity']",
    "css" : {
      "width" : 110.0
    }
  }, {
    "selector" : "node[bel_function_type = 'cellSurfaceExpression']",
    "css" : {
      "width" : 110.0
    }
  }, {
    "selector" : "node[bel_function_type = 'chaperoneActivity']",
    "css" : {
      "width" : 110.0
    }
  }, {
    "selector" : "node[bel_function_type = 'transcriptionalActivity']",
    "css" : {
      "width" : 110.0
    }
  }, {
    "selector" : "node[bel_function_type = 'phosphataseActivity']",
    "css" : {
      "width" : 110.0
    }
  }, {
    "selector" : "node[bel_function_type = 'complexAbundance']",
    "css" : {
      "width" : 150.0
    }
  }, {
    "selector" : "node[linked]",
    "css" : {
      "border-opacity" : 1.0
    }
  }, {
    "selector" : "node[!linked]",
    "css" : {
      "border-opacity" : 0.39215686274509803
    }
  }, {
    "selector" : "node[linked]",
    "css" : {
      "text-opacity" : 1.0
    }
  }, {
    "selector" : "node[!linked]",
    "css" : {
      "text-opacity" : 0.39215686274509803
    }
  }, {
    "selector" : "node[linked]",
    "css" : {
      "background-opacity" : 1.0
    }
  }, {
    "selector" : "node[!linked]",
    "css" : {
      "background-opacity" : 0.39215686274509803
    }
  }, {
    "selector" : "node[bel_function_type = 'transcriptionalActivity']",
    "css" : {
      "width" : 100.0,
      "height" : 100.0
    }
  }, {
    "selector" : "node[bel_function_type = 'transportActivity']",
    "css" : {
      "shape" : "roundrectangle"
    }
  }, {
    "selector" : "node[bel_function_type = 'cellSecretion']",
    "css" : {
      "shape" : "v"
    }
  }, {
    "selector" : "node[bel_function_type = 'catalyticActivity']",
    "css" : {
      "shape" : "hexagon"
    }
  }, {
    "selector" : "node[bel_function_type = 'reaction']",
    "css" : {
      "shape" : "parallelogram"
    }
  }, {
    "selector" : "node[bel_function_type = 'pathology']",
    "css" : {
      "shape" : "parallelogram"
    }
  }, {
    "selector" : "node[bel_function_type = 'gtpBoundActivity']",
    "css" : {
      "shape" : "hexagon"
    }
  }, {
    "selector" : "node[bel_function_type = 'translocation']",
    "css" : {
      "shape" : "v"
    }
  }, {
    "selector" : "node[bel_function_type = 'kinaseActivity']",
    "css" : {
      "shape" : "hexagon"
    }
  }, {
    "selector" : "node[bel_function_type = 'proteinAbundance']",
    "css" : {
      "shape" : "roundrectangle"
    }
  }, {
    "selector" : "node[bel_function_type = 'peptidaseActivity']",
    "css" : {
      "shape" : "hexagon"
    }
  }, {
    "selector" : "node[bel_function_type = 'abundance']",
    "css" : {
      "shape" : "diamond"
    }
  }, {
    "selector" : "node[bel_function_type = 'degradation']",
    "css" : {
      "shape" : "roundrectangle"
    }
  }, {
    "selector" : "node[bel_function_type = 'compositeAbundance']",
    "css" : {
      "shape" : "roundrectangle"
    }
  }, {
    "selector" : "node[bel_function_type = 'molecularActivity']",
    "css" : {
      "shape" : "hexagon"
    }
  }, {
    "selector" : "node[bel_function_type = 'rnaAbundance']",
    "css" : {
      "shape" : "ellipse"
    }
  }, {
    "selector" : "node[bel_function_type = 'biologicalProcess']",
    "css" : {
      "shape" : "parallelogram"
    }
  }, {
    "selector" : "node[bel_function_type = 'microRNAAbundance']",
    "css" : {
      "shape" : "ellipse"
    }
  }, {
    "selector" : "node[bel_function_type = 'geneAbundance']",
    "css" : {
      "shape" : "ellipse"
    }
  }, {
    "selector" : "node[bel_function_type = 'ribosylationActivity']",
    "css" : {
      "shape" : "hexagon"
    }
  }, {
    "selector" : "node[bel_function_type = 'cellSurfaceExpression']",
    "css" : {
      "shape" : "v"
    }
  }, {
    "selector" : "node[bel_function_type = 'chaperoneActivity']",
    "css" : {
      "shape" : "hexagon"
    }
  }, {
    "selector" : "node[bel_function_type = 'transcriptionalActivity']",
    "css" : {
      "shape" : "hexagon"
    }
  }, {
    "selector" : "node[bel_function_type = 'phosphataseActivity']",
    "css" : {
      "shape" : "hexagon"
    }
  }, {
    "selector" : "node[bel_function_type = 'complexAbundance']",
    "css" : {
      "shape" : "roundrectangle"
    }
  }, {
    "selector" : "node[ id = '107' ]",
    "css" : { }
  }, {
    "selector" : "node:selected",
    "css" : {
      "background-color" : "rgb(255,255,0)"
    }
  }, {
    "selector" : "edge",
    "css" : {
      "font-size" : 10,
      "target-arrow-color" : "rgb(0,0,0)",
      "source-arrow-shape" : "none",
      "width" : 1.0,
      "content" : "",
      "opacity" : 1.0,
      "color" : "rgb(0,0,0)",
      "line-color" : "rgb(0,0,0)",
      "line-style" : "solid",
      "target-arrow-shape" : "none",
      "text-opacity" : 1.0,
      "source-arrow-color" : "rgb(0,0,0)",
      "font-family" : "Dialog",
      "font-weight" : "normal"
    }
  }, {
    "selector" : "edge[interaction = 'transcribedTo']",
    "css" : {
      "line-style" : "dashed"
    }
  }, {
    "selector" : "edge[interaction = 'causesNoChange']",
    "css" : {
      "line-style" : "dashed"
    }
  }, {
    "selector" : "edge[interaction = 'translatedTo']",
    "css" : {
      "line-style" : "dashed"
    }
  }, {
    "selector" : "edge[interaction = 'hasProduct']",
    "css" : {
      "line-style" : "dashed"
    }
  }, {
    "selector" : "edge[interaction = 'positiveCorrelation']",
    "css" : {
      "line-style" : "dotted"
    }
  }, {
    "selector" : "edge[interaction = 'reactantIn']",
    "css" : {
      "line-style" : "dashed"
    }
  }, {
    "selector" : "edge[interaction = 'subProcessOf']",
    "css" : {
      "line-style" : "dashed"
    }
  }, {
    "selector" : "edge[interaction = 'hasMember']",
    "css" : {
      "line-style" : "dashed"
    }
  }, {
    "selector" : "edge[interaction = 'orthologous']",
    "css" : {
      "line-style" : "dashed"
    }
  }, {
    "selector" : "edge[interaction = 'association']",
    "css" : {
      "line-style" : "dotted"
    }
  }, {
    "selector" : "edge[interaction = 'actsIn']",
    "css" : {
      "line-style" : "dashed"
    }
  }, {
    "selector" : "edge[interaction = 'includes']",
    "css" : {
      "line-style" : "dashed"
    }
  }, {
    "selector" : "edge[interaction = 'hasVariant']",
    "css" : {
      "line-style" : "dashed"
    }
  }, {
    "selector" : "edge[interaction = 'directlyDecreases']",
    "css" : {
      "line-style" : "solid"
    }
  }, {
    "selector" : "edge[interaction = 'negativeCorrelation']",
    "css" : {
      "line-style" : "dotted"
    }
  }, {
    "selector" : "edge[interaction = 'directlyIncreases']",
    "css" : {
      "line-style" : "solid"
    }
  }, {
    "selector" : "edge[interaction = 'prognosticBiomarkerFor']",
    "css" : {
      "line-style" : "dotted"
    }
  }, {
    "selector" : "edge[interaction = 'hasComponent']",
    "css" : {
      "line-style" : "dashed"
    }
  }, {
    "selector" : "edge[interaction = 'hasModification']",
    "css" : {
      "line-style" : "dashed"
    }
  }, {
    "selector" : "edge[interaction = 'subProcess_OF']",
    "css" : {
      "line-style" : "dashed"
    }
  }, {
    "selector" : "edge[interaction = 'isA']",
    "css" : {
      "line-style" : "dashed"
    }
  }, {
    "selector" : "edge[interaction = 'rateLimitingStepOf']",
    "css" : {
      "line-style" : "dashed"
    }
  }, {
    "selector" : "edge[interaction = 'biomarkerFor']",
    "css" : {
      "line-style" : "dotted"
    }
  }, {
    "selector" : "edge[interaction = 'translocates']",
    "css" : {
      "line-style" : "dashed"
    }
  }, {
    "selector" : "edge[interaction = 'translatedTo']",
    "css" : {
      "target-arrow-shape" : "diamond"
    }
  }, {
    "selector" : "edge[interaction = 'hasProduct']",
    "css" : {
      "target-arrow-shape" : "diamond"
    }
  }, {
    "selector" : "edge[interaction = 'reactantIn']",
    "css" : {
      "target-arrow-shape" : "diamond"
    }
  }, {
    "selector" : "edge[interaction = 'orthologous']",
    "css" : {
      "target-arrow-shape" : "none"
    }
  }, {
    "selector" : "edge[interaction = 'actsIn']",
    "css" : {
      "target-arrow-shape" : "diamond"
    }
  }, {
    "selector" : "edge[interaction = 'association']",
    "css" : {
      "target-arrow-shape" : "none"
    }
  }, {
    "selector" : "edge[interaction = 'hasVariant']",
    "css" : {
      "target-arrow-shape" : "diamond"
    }
  }, {
    "selector" : "edge[interaction = 'directlyIncreases']",
    "css" : {
      "target-arrow-shape" : "triangle"
    }
  }, {
    "selector" : "edge[interaction = 'unknown']",
    "css" : {
      "target-arrow-shape" : "triangle"
    }
  }, {
    "selector" : "edge[interaction = 'increases']",
    "css" : {
      "target-arrow-shape" : "triangle"
    }
  }, {
    "selector" : "edge[interaction = 'prognosticBiomarkerFor']",
    "css" : {
      "target-arrow-shape" : "triangle"
    }
  }, {
    "selector" : "edge[interaction = 'hasComponent']",
    "css" : {
      "target-arrow-shape" : "diamond"
    }
  }, {
    "selector" : "edge[interaction = 'hasModification']",
    "css" : {
      "target-arrow-shape" : "diamond"
    }
  }, {
    "selector" : "edge[interaction = 'subProcess_OF']",
    "css" : {
      "target-arrow-shape" : "diamond"
    }
  }, {
    "selector" : "edge[interaction = 'isA']",
    "css" : {
      "target-arrow-shape" : "diamond"
    }
  }, {
    "selector" : "edge[interaction = 'rateLimitingStepOf']",
    "css" : {
      "target-arrow-shape" : "triangle"
    }
  }, {
    "selector" : "edge[interaction = 'biomarkerFor']",
    "css" : {
      "target-arrow-shape" : "triangle"
    }
  }, {
    "selector" : "edge[interaction = 'translocates']",
    "css" : {
      "target-arrow-shape" : "diamond"
    }
  }, {
    "selector" : "edge[interaction = 'analogous']",
    "css" : {
      "target-arrow-shape" : "triangle"
    }
  }, {
    "selector" : "edge[interaction = 'transcribedTo']",
    "css" : {
      "target-arrow-shape" : "diamond"
    }
  }, {
    "selector" : "edge[interaction = 'positiveCorrelation']",
    "css" : {
      "target-arrow-shape" : "triangle"
    }
  }, {
    "selector" : "edge[interaction = 'subProcessOf']",
    "css" : {
      "target-arrow-shape" : "diamond"
    }
  }, {
    "selector" : "edge[interaction = 'hasMember']",
    "css" : {
      "target-arrow-shape" : "diamond"
    }
  }, {
    "selector" : "edge[interaction = 'includes']",
    "css" : {
      "target-arrow-shape" : "diamond"
    }
  }, {
    "selector" : "edge[interaction = 'directlyDecreases']",
    "css" : {
      "target-arrow-shape" : "circle"
    }
  }, {
    "selector" : "edge[interaction = 'negativeCorrelation']",
    "css" : {
      "target-arrow-shape" : "circle"
    }
  }, {
    "selector" : "edge[interaction = 'decreases']",
    "css" : {
      "target-arrow-shape" : "circle"
    }
  }, {
    "selector" : "edge[interaction = 'positiveCorrelation']",
    "css" : {
      "source-arrow-shape" : "triangle"
    }
  }, {
    "selector" : "edge[interaction = 'hasMember']",
    "css" : {
      "source-arrow-shape" : "none"
    }
  }, {
    "selector" : "edge[interaction = 'orthologous']",
    "css" : {
      "source-arrow-shape" : "none"
    }
  }, {
    "selector" : "edge[interaction = 'association']",
    "css" : {
      "source-arrow-shape" : "none"
    }
  }, {
    "selector" : "edge[interaction = 'negativeCorrelation']",
    "css" : {
      "source-arrow-shape" : "circle"
    }
  }, {
    "selector" : "edge[interaction = 'translatedTo']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'hasProduct']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'reactantIn']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'orthologous']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'actsIn']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'association']",
    "css" : {
      "width" : 1.0
    }
  }, {
    "selector" : "edge[interaction = 'hasVariant']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'directlyIncreases']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'increases']",
    "css" : {
      "width" : 1.0
    }
  }, {
    "selector" : "edge[interaction = 'prognosticBiomarkerFor']",
    "css" : {
      "width" : 1.0
    }
  }, {
    "selector" : "edge[interaction = 'hasComponent']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'hasModification']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'subProcess_OF']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'isA']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'rateLimitingStepOf']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'biomarkerFor']",
    "css" : {
      "width" : 1.0
    }
  }, {
    "selector" : "edge[interaction = 'translocates']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'transcribedTo']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'causesNoChange']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'positiveCorrelation']",
    "css" : {
      "width" : 1.0
    }
  }, {
    "selector" : "edge[interaction = 'subProcessOf']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'hasMember']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'includes']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'directlyDecreases']",
    "css" : {
      "width" : 3.0
    }
  }, {
    "selector" : "edge[interaction = 'negativeCorrelation']",
    "css" : {
      "width" : 1.0
    }
  }, {
    "selector" : "edge[interaction = 'decreases']",
    "css" : {
      "width" : 1.0
    }
  }, {
    "selector" : "edge[linked]",
    "css" : {
      "opacity" : 1.0
    }
  }, {
    "selector" : "edge[!linked]",
    "css" : {
      "opacity" : 0.19607843137254902
    }
  }, {
    "selector" : "edge[interaction = 'transcribedTo']",
    "css" : {
      "line-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'translatedTo']",
    "css" : {
      "line-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'hasProduct']",
    "css" : {
      "line-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'reactantIn']",
    "css" : {
      "line-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'subProcessOf']",
    "css" : {
      "line-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'hasMember']",
    "css" : {
      "line-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'orthologous']",
    "css" : {
      "line-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'actsIn']",
    "css" : {
      "line-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'includes']",
    "css" : {
      "line-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'hasVariant']",
    "css" : {
      "line-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'hasComponent']",
    "css" : {
      "line-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'hasModification']",
    "css" : {
      "line-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'subProcess_OF']",
    "css" : {
      "line-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'isA']",
    "css" : {
      "line-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'translocates']",
    "css" : {
      "line-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'transcribedTo']",
    "css" : {
      "line-color" : "rgb(153,153,153)",
      "target-arrow-color" : "rgb(153,153,153)",
      "source-arrow-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'translatedTo']",
    "css" : {
      "line-color" : "rgb(153,153,153)",
      "target-arrow-color" : "rgb(153,153,153)",
      "source-arrow-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'hasProduct']",
    "css" : {
      "line-color" : "rgb(153,153,153)",
      "target-arrow-color" : "rgb(153,153,153)",
      "source-arrow-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'reactantIn']",
    "css" : {
      "line-color" : "rgb(153,153,153)",
      "target-arrow-color" : "rgb(153,153,153)",
      "source-arrow-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'subProcessOf']",
    "css" : {
      "line-color" : "rgb(153,153,153)",
      "target-arrow-color" : "rgb(153,153,153)",
      "source-arrow-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'hasMember']",
    "css" : {
      "line-color" : "rgb(153,153,153)",
      "target-arrow-color" : "rgb(153,153,153)",
      "source-arrow-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'orthologous']",
    "css" : {
      "line-color" : "rgb(153,153,153)",
      "target-arrow-color" : "rgb(153,153,153)",
      "source-arrow-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'actsIn']",
    "css" : {
      "line-color" : "rgb(153,153,153)",
      "target-arrow-color" : "rgb(153,153,153)",
      "source-arrow-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'includes']",
    "css" : {
      "line-color" : "rgb(153,153,153)",
      "target-arrow-color" : "rgb(153,153,153)",
      "source-arrow-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'hasVariant']",
    "css" : {
      "line-color" : "rgb(153,153,153)",
      "target-arrow-color" : "rgb(153,153,153)",
      "source-arrow-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'hasComponent']",
    "css" : {
      "line-color" : "rgb(153,153,153)",
      "target-arrow-color" : "rgb(153,153,153)",
      "source-arrow-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'hasModification']",
    "css" : {
      "line-color" : "rgb(153,153,153)",
      "target-arrow-color" : "rgb(153,153,153)",
      "source-arrow-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'subProcess_OF']",
    "css" : {
      "line-color" : "rgb(153,153,153)",
      "target-arrow-color" : "rgb(153,153,153)",
      "source-arrow-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'isA']",
    "css" : {
      "line-color" : "rgb(153,153,153)",
      "target-arrow-color" : "rgb(153,153,153)",
      "source-arrow-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge[interaction = 'translocates']",
    "css" : {
      "line-color" : "rgb(153,153,153)",
      "target-arrow-color" : "rgb(153,153,153)",
      "source-arrow-color" : "rgb(153,153,153)"
    }
  }, {
    "selector" : "edge:selected",
    "css" : {
      "line-color" : "rgb(255,0,0)"
    }
  } ]
}, {
  "format_version" : "1.0",
  "generated_by" : "cytoscape-3.7.2",
  "target_cytoscapejs_version" : "~2.1",
  "title" : "Marquee",
  "style" : [ {
    "selector" : "node",
    "css" : {
      "border-width" : 10.0,
      "height" : 20.0,
      "shape" : "ellipse",
      "text-valign" : "bottom",
      "text-halign" : "center",
      "border-color" : "rgb(255,255,255)",
      "color" : "rgb(102,102,102)",
      "font-size" : 12,
      "font-family" : "SansSerif",
      "font-weight" : "normal",
      "width" : 20.0,
      "border-opacity" : 1.0,
      "text-opacity" : 1.0,
      "background-color" : "rgb(0,204,255)",
      "background-opacity" : 1.0,
      "content" : "data(name)"
    }
  }, {
    "selector" : "node[ id = '107' ]",
    "css" : { }
  }, {
    "selector" : "node:selected",
    "css" : {
      "background-color" : "rgb(255,0,102)"
    }
  }, {
    "selector" : "edge",
    "css" : {
      "font-size" : 8,
      "target-arrow-color" : "rgb(255,255,255)",
      "source-arrow-shape" : "none",
      "width" : 2.0,
      "opacity" : 1.0,
      "color" : "rgb(102,102,102)",
      "line-color" : "rgb(255,255,255)",
      "line-style" : "dashed",
      "target-arrow-shape" : "triangle",
      "text-opacity" : 1.0,
      "source-arrow-color" : "rgb(255,255,255)",
      "font-family" : "SansSerif",
      "font-weight" : "normal",
      "content" : "data(interaction)"
    }
  }, {
    "selector" : "edge:selected",
    "css" : {
      "line-color" : "rgb(255,0,0)"
    }
  } ]
}, {
  "format_version" : "1.0",
  "generated_by" : "cytoscape-3.7.2",
  "target_cytoscapejs_version" : "~2.1",
  "title" : "Sample2",
  "style" : [ {
    "selector" : "node",
    "css" : {
      "border-width" : 15.0,
      "height" : 50.0,
      "shape" : "ellipse",
      "text-valign" : "center",
      "text-halign" : "right",
      "border-color" : "rgb(255,255,255)",
      "color" : "rgb(102,102,102)",
      "font-size" : 20,
      "font-family" : "HelveticaNeue-Light",
      "font-weight" : "normal",
      "width" : 50.0,
      "border-opacity" : 1.0,
      "text-opacity" : 1.0,
      "background-color" : "rgb(58,127,182)",
      "background-opacity" : 1.0,
      "content" : "data(name)"
    }
  }, {
    "selector" : "node[ id = '107' ]",
    "css" : { }
  }, {
    "selector" : "node:selected",
    "css" : {
      "background-color" : "rgb(255,255,0)"
    }
  }, {
    "selector" : "edge",
    "css" : {
      "font-size" : 10,
      "target-arrow-color" : "rgb(0,0,0)",
      "source-arrow-shape" : "none",
      "width" : 20.0,
      "content" : "",
      "opacity" : 1.0,
      "color" : "rgb(0,0,0)",
      "line-color" : "rgb(255,255,255)",
      "line-style" : "solid",
      "target-arrow-shape" : "none",
      "text-opacity" : 1.0,
      "source-arrow-color" : "rgb(0,0,0)",
      "font-family" : "SansSerif",
      "font-weight" : "normal"
    }
  }, {
    "selector" : "edge:selected",
    "css" : {
      "line-color" : "rgb(255,0,0)"
    }
  } ]
}, {
  "format_version" : "1.0",
  "generated_by" : "cytoscape-3.7.2",
  "target_cytoscapejs_version" : "~2.1",
  "title" : "Directed",
  "style" : [ {
    "selector" : "node",
    "css" : {
      "border-width" : 5.0,
      "height" : 45.0,
      "shape" : "ellipse",
      "text-valign" : "center",
      "text-halign" : "center",
      "border-color" : "rgb(145,145,145)",
      "color" : "rgb(51,153,255)",
      "font-size" : 8,
      "font-family" : "SansSerif",
      "font-weight" : "normal",
      "width" : 45.0,
      "border-opacity" : 1.0,
      "text-opacity" : 1.0,
      "background-color" : "rgb(255,255,255)",
      "background-opacity" : 1.0,
      "content" : "data(name)"
    }
  }, {
    "selector" : "node[ id = '107' ]",
    "css" : { }
  }, {
    "selector" : "node:selected",
    "css" : {
      "background-color" : "rgb(255,0,102)"
    }
  }, {
    "selector" : "edge",
    "css" : {
      "font-size" : 12,
      "target-arrow-color" : "rgb(204,204,204)",
      "source-arrow-shape" : "none",
      "width" : 5.0,
      "opacity" : 1.0,
      "color" : "rgb(51,153,255)",
      "line-color" : "rgb(204,204,204)",
      "line-style" : "solid",
      "target-arrow-shape" : "triangle",
      "text-opacity" : 1.0,
      "source-arrow-color" : "rgb(204,204,204)",
      "font-family" : "SansSerif",
      "font-weight" : "normal",
      "content" : "data(interaction)"
    }
  }, {
    "selector" : "edge:selected",
    "css" : {
      "line-color" : "rgb(255,0,0)"
    }
  } ]
}, {
  "format_version" : "1.0",
  "generated_by" : "cytoscape-3.7.2",
  "target_cytoscapejs_version" : "~2.1",
  "title" : "Sample1",
  "style" : [ {
    "selector" : "node",
    "css" : {
      "border-width" : 0.0,
      "height" : 25.0,
      "shape" : "ellipse",
      "text-valign" : "center",
      "text-halign" : "center",
      "border-color" : "rgb(0,0,0)",
      "color" : "rgb(51,51,51)",
      "font-size" : 10,
      "font-family" : "Dialog",
      "font-weight" : "normal",
      "width" : 25.0,
      "border-opacity" : 1.0,
      "text-opacity" : 1.0,
      "background-color" : "rgb(127,205,187)",
      "background-opacity" : 1.0,
      "content" : "data(name)"
    }
  }, {
    "selector" : "node[ id = '107' ]",
    "css" : { }
  }, {
    "selector" : "node:selected",
    "css" : {
      "background-color" : "rgb(255,255,0)"
    }
  }, {
    "selector" : "edge",
    "css" : {
      "font-size" : 10,
      "target-arrow-color" : "rgb(0,0,0)",
      "source-arrow-shape" : "none",
      "width" : 1.0,
      "opacity" : 1.0,
      "color" : "rgb(51,51,51)",
      "line-color" : "rgb(153,153,153)",
      "line-style" : "solid",
      "target-arrow-shape" : "none",
      "text-opacity" : 1.0,
      "source-arrow-color" : "rgb(0,0,0)",
      "font-family" : "Dialog",
      "font-weight" : "normal",
      "content" : "data(interaction)"
    }
  }, {
    "selector" : "edge[interaction = 'pp']",
    "css" : {
      "line-style" : "solid"
    }
  }, {
    "selector" : "edge[interaction = 'pd']",
    "css" : {
      "line-style" : "dashed"
    }
  }, {
    "selector" : "edge:selected",
    "css" : {
      "line-color" : "rgb(255,0,0)"
    }
  } ]
} ]