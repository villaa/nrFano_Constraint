function toggle_visibility(id) {
   var e = document.getElementById(id);
   if(e.style.display == 'block')
      e.style.display = 'none';
   else
      e.style.display = 'block';
}
function remove_visibility(id) {
   var e = document.getElementById(id);
   e.style.display = 'none';
}
function restore_visibility(id) {
   var e = document.getElementById(id);
   e.style.display = 'block';
}
