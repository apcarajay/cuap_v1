$(document).ready(function() {

	// Funtion to get checked radio's value.
	$('#radio_value').click(function() {
		$('#result').empty();
		var value = $("form input[type='radio']:checked").val();
		var value1 = $('#radio_value').val();
		if($("form input[type='radio']").is(':checked')) {
			//$('#result').append(value);
			Shiny.onInputChange("radioval",value);
			//alert(value1);
			//Shiny.onInputChange("radio_value",value1);
		}else{
			//value = "zero";
			//Shiny.onInputChange("radioval",value);
			alert(" Please click any option! ");
			location.reload();
		}

	});

	// Get value Onchange radio function.
	$('input:radio').change(function(){
		var value = $("form input[type='radio']:checked").val();
		Shiny.onInputChange("radioval",value);
		//alert("Value of Changed Radio is : " +value);
	});


	// file upload
	$('#inputfile').change(function(event) {
        var filename = $('#inputfile').val();
        if (filename.substring(3,11) == 'fakepath') {
            filename = filename.substring(12);
        } // Remove c:\fake at beginning from localhost chrome
       
       	var accept = $(event.target).attr('accept');
	    //<input name='fileImage' accept='.png,.jpg,.jpeg' type='file'  class='fileImage' />
	    value = $(event.target).val(), ext = value.substring(value.lastIndexOf('.'));
	    //alert(accept); // .fasta,.txt
	    //alert(value);  // C:fakepath/sample.fasta
	    //alert(ext); // .fasta
	    //alert(accept.indexOf(ext));
	    if(ext && accept.indexOf(ext) == -1){
	    	alert("INVALID FILE TYPE. ONLY ACCEPTS .FASTA OR .TXT FILE");
	        location.reload();

	    }else{
	    	//alert("FILENAME: "+filename);
	    	 //var myFile = $('#inputfile').prop('files')[0];
        	//alert(myFile);
        	  var file = event.target.files[0];
			  var reader = new FileReader();
			  reader.onload = function(event) {
			    // The file's text will be printed here
			    var data = event.target.result;
			    var lines = data.split("\n");
			    //alert(data);
			    //alert(lines.length);
			    //lines=lines.replace(/\t/g, '');
			    //alert(lines);
			    var inv=0;
			    var index=[];
			    for (var i = 0, len = lines.length; i < len; i++) {
		            //alert(lines[i]);
		            var perline = lines[i].split("");

		            for(var j = 0, jlen = perline.length; j < jlen; j++){
		            	//alert(" - "+ perline[j]);
		            	if(perline[j] == ">" && j==0){
		            		var inv = -1;
		            		//alert("Index:" + index.push(inv) + lines[i] + i);
		            		break;
		            	}else{
		            		var inv = 0;
		            	}
		            }
		            if(inv== -1){
		            	index.push(-1);
		            	//alert("Index:" + index.push(-1) + lines[i] + i);
		            }else{
		            	index.push(0);
		            	//alert("Index:" + index.push(0)  + lines[i] + i);
		            }
		        }
		        //alert(index);
		        //alert(index[0]);
		        //alert(index.length);

		        count = 0;
		        for (var k = 0, klen = lines.length; k < klen; k++) {
		            //alert(lines[i]);
		            
		            //alert("index["+k+"] is "+ index[k]);
		            if(index[k] == -1){
		            	//alert("dont read");
		            	//break;
		            //}else if(k > index.length && eachline[0] == ">"){
		            	//alert("beyond index");

		            }else{
		            	var eachline = lines[k].split("");
		            	n = k + 1;
		            	//eachline = eachlines[i].split(' ').join('');
		            	//var eachline = eachline1.replace(/\s+/g, ' ').trim();

		            	for(var m = 0, mlen = eachline.length; m < mlen; m++){
			            	//alert(" - "+ perline[j]);

			            	if(eachline[m] == "A" || eachline[m] == "T" || eachline[m] == "G" || eachline[m] == "C" || eachline[m] == "a" || eachline[m] == "t" || eachline[m] == "g" || eachline[m] == "c"){
			            					            		
			            	}else if(eachline[m] == " " || eachline[m] == "\t" || eachline[m] == "\n" || eachline[m] == "\r"){
			            		
			            	}else{
			            		alert("Error: Input contains invalid base " + eachline[m] + " at line "+ n);
			            		//break;
			            		location.reload(true);
			            		exit();
			            		//return false;
			            	}
			            }
		            }
		            
		        }

			  };

			  reader.readAsText(file);

	    }
   });



	// Funtion to reset or clear selection.
	$('#form_reset').click(function() {
		//alert("Form reset");
		$('#result').empty();
		$("input:radio").attr("checked", false);
		location.reload();
	});

	// Funtion to back input selection.
	$('#backtoinput').click(function() {
		//alert("Form reset");
		$('#result').empty();
		$("input:radio").attr("checked", false);
		location.reload();
	});


	$('#inputfile').filestyle({ 
		buttonName : 'btn-success',		 
		buttonText : ' Open'	 
	});

	

	$('#myTable').DataTable( {
    responsive: true
} );
});