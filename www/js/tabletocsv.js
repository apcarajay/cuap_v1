
$(document).ready(function() {
	$('#btnDownload').click(function() {
		$("#mytable").tableToCSV({
			filename: 'Summaryofdata'
		});
	});
});
