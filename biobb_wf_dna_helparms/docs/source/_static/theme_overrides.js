$(document).ready(function() {
	var path = '';
	
	if($('.wy-side-nav-search .icon.icon-home').attr('href') == '../index.html') path = '../';
	else if($('.wy-side-nav-search .icon.icon-home').attr('href') == '../../index.html') path = '../../';
		
	$('.wy-side-nav-search .icon.icon-home').html('<img src="' + path + '_static/logo.png" class="logo" alt="Logo">');

	$('table.dataframe').parent().css('overflow-x', 'auto');
	$('table.dataframe').parent().css('margin-bottom', '20px');
	$('table.dataframe').css('border-color', '#999');
	$('table.dataframe').css('color', '#555');
	$('table.dataframe td, table.dataframe th').css('padding', '5px 10px');
});