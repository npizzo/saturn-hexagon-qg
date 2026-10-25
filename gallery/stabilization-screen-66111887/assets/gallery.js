
document.querySelectorAll('[data-frame-viewer]').forEach(function(viewer){
  const data=JSON.parse(viewer.querySelector('script[type="application/json"]').textContent);
  const slider=viewer.querySelector('input[type="range"]');
  const image=viewer.querySelector('img');
  const label=viewer.querySelector('[data-frame-label]');
  function update(){const i=Number(slider.value), f=data[i];image.src=f.image;image.alt='Relative vorticity and speed at '+f.time_s+' model seconds';label.textContent='Saved frame '+(i+1)+' / '+data.length+' · '+(f.time_s/86400).toFixed(3)+' days · '+f.time_s.toLocaleString()+' s';}
  slider.addEventListener('input',update);
  viewer.querySelector('[data-previous]').addEventListener('click',function(){slider.value=Math.max(0,Number(slider.value)-1);update();});
  viewer.querySelector('[data-next]').addEventListener('click',function(){slider.value=Math.min(data.length-1,Number(slider.value)+1);update();});
  update();
});
