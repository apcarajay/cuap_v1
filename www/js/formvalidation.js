// for abled/disabled radio button 
<script>

	window.onload=function() {

    // first, disable all the input fields
    document.forms[0].elements["sequence"].disabled=true;
    document.forms[0].elements["accessno"].disabled=true;
    document.forms[0].elements["inputfile"].disabled=true;

    // next, attach the click event handler to the radio buttons
    var radios = document.forms[0].elements["group1"];
	    for (var i = [0]; i < radios.length; i++){
	  		radios[i].onclick=radioClicked;
    	}
		    
	}

		function radioClicked() {

		  // find out which radio button was clicked and
		  // disable/enable appropriate input elements
		  switch(this.value) {
		    case "one" :
		       //document.location.reload(true);
		       document.forms[0].elements["sequence"].disabled=false;
		       document.forms[0].elements["accessno"].disabled=true;
		       document.forms[0].elements["inputfile"].disabled=true;
			   break;
		    case "two" :
		       //document.location.reload(true);
		       document.forms[0].elements["accessno"].disabled=false;
		       document.forms[0].elements["sequence"].disabled=true;
		       document.forms[0].elements["inputfile"].disabled=true;
		       break;
		    case "three" :
		       //document.location.reload(true);
		       document.forms[0].elements["inputfile"].disabled=false;
		       document.forms[0].elements["sequence"].disabled=true;
		       document.forms[0].elements["accessno"].disabled=true;
		       break;
		  }
		}


		var valid = false;

		function validate_fileupload(input_element)
		{
		    var el = document.getElementById("feedback");
		    var fileName = input_element.value;
		    var allowed_extensions = new Array("fasta","txt");
		    var file_extension = fileName.split('.').pop(); 
		    for(var i = 0; i < allowed_extensions.length; i++)
		    {
		        if(allowed_extensions[i]==file_extension)
		        {
		            valid = true; // valid file extension
		            el.innerHTML = "";
		            return;
		        }
		    }
		    el.innerHTML="Invalid file";
		    el.form.reset();
            el.focus();
		    valid = false;

		}

		function ValidateInputForm(){
			var sequence = document.getElementById("sequence");           
			var accessno = document.getElementById("accessno");
			var one = document.getElementById("one").checked;
			var two = document.getElementById("two").checked;
			var three = document.getElementById("three").checked;
			var lines = sequence.value.split("\n");
			var countseq = 0;
			//var pattern = ^[ATGCatgc]+$;
			//var res = pattern.test(sequence.value);


    		if (sequence.value == "" && one == true){
		        window.alert("Please enter your sequence(s).");
		        sequence.focus();
		        document.location.reload(true);
		        return false;
		    }

		    if (sequence.value != "" && one == true){
		        //window.alert(sequence.value);
		        for(var j=0; j < lines.length; j++){
		        	//window.alert("Line"+ j + "is : " + lines[j]);
		        	window.alert(lines[j].indexOf(">"))
		        	if (lines[j].indexOf(">",0) < 0){
		        		window.alert(lines[j] + ";  check if valid dna(ATGC");
		        		var res = pattern.test(lines[j]);
		        		window.alert(res);
		        	}else{
		        		 countseq++;
		        		 window.alert("count is "+ countseq);
		        	}

		        	if(countseq > 0){
		        		return true;
		        	}
		        }
		        sequence.focus();
		        document.location.reload(true);
		        return false;
		    }

		    if (accessno.value == "" && two == true){
		        window.alert("Please enter your accession number(s).");
		        accessno.focus();
		        document.location.reload(true);
		        return false;
		    }   
  			return true;
		}
		
		function resetForm(){
			//document.getElementById("picker").reset();
			document.location.reload(true);
		}
		
</script>