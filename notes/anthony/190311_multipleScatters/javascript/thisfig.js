var toggleShow = function(id_fig,id_str) {
  var svgDoc = document.getElementById(id_fig).contentDocument;
  var thisLine = svgDoc.getElementById(id_str);
  console.log(thisLine);
  thisLine.classList.toggle('hidden');
};
