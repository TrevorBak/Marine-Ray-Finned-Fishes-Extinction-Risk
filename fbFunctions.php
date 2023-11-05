<?php	
	# custom functions for extracting information from fishbase
	# Written by Noel A. Heim, Tufts University, noel.heim@tufts.edu
	# Last modified 20 November 2020
	
	// extract reference from string and return ref id
	function getRef($fbstring) {
		$x = array();
		$y = preg_match("/(\(Ref. )(.*?)\)/s", $fbstring, $x);
		return(end($x));
	}
	
	// extract reference from string and return ref id
	function getEcoUrl($pagestring) {
		$eco_base = "http://www.fishbase.us/Ecology/FishEcologySummary.php";
		$x = array();
		$y = preg_match("/(FishEcologySummary\.php\?)(.*?)( alt\=\'Ecology for)/s", $pagestring, $x);
		return($eco_base."?".trim(rtrim($x[2],"\'")));
	}
	
	// function to extract phylogenetic diversity
	function getPhyloDiv($pagestring, $taxonId) {
		$myReturn = array('pdiv' => NULL, 'pdivref' => NULL, 'tlevel' => NULL, 'tlevelse' => NULL, 'tlevelref' => NULL, 
			'resil' => NULL, 'resil_ref' => NULL, 'resil_K' => NULL, 'resil_tm' => NULL, 'resil_tmax' => NULL, 'resil_fec' => NULL,
			'vulnerab' => NULL, 'vulnerab_ref' => NULL, 'price' => NULL, 'price_ref' => NULL,
			'main_ref' => NULL);
		
		//PRICE
		$beginstr = "<!-- Start price-->";
		$endstr = "<\/div>";
		$x = array();
		$y = preg_match("/(".$beginstr.")(.*?)(".$endstr.")/s", $pagestring, $x);

		if(count($x)==4) {
			$selectedstr = $x[2];
			// get value
			if(preg_match("/unknown<\/a>/i", $selectedstr)==1) $myReturn['price'] = 'unknown';
			else if(preg_match("/low<\/a>/i", $selectedstr)==1) $myReturn['price'] = 'low';
			else if(preg_match("/medium<\/a>/i", $selectedstr)==1) $myReturn['price'] = 'medium';
			else if(preg_match("/very high<\/a>/i", $selectedstr)==1) $myReturn['price'] = 'very high';
			else if(preg_match("/high<\/a>/i", $selectedstr)==1) $myReturn['price'] = 'high';
			
			// get reference
			$beginstr = "references/FBRefSummary.php\?ID=[0-9]+'>";
			$endstr = "</a>";
			$ref = array();
			$y = preg_match("#(".$beginstr.")(.*?)(".$endstr.")#i", $selectedstr, $ref);			
			if(count($ref) == 4) $myReturn['price_ref'] = trim($ref[2]);
			unset($ref);
		} 
		
		
		//MAIN REFERENCE
		$beginstr = "<!-- Main Reference -->";
		$endstr = "<\/p>";
		$x = array();
		$y = preg_match("/(".$beginstr.")(.*?)(".$endstr.")/s", $pagestring, $x);

		if(count($x)==4) {
			$selectedstr = $x[2]; 
			// get reference
			$beginstr = "references/FBRefSummary.php\?ID=[0-9]+'>";
			$endstr = "</a>";
			$ref = array();
			$y = preg_match("#(".$beginstr.")(.*?)(".$endstr.")#i", $selectedstr, $ref);			
			if(count($ref) == 4) $myReturn['main_ref'] = trim($ref[2]);
			unset($ref);
		} 
		
		//VULNERABILITY
		$beginstr = "<!-- Start vulnerability-->";
		$endstr = "<\/div>";
		$x = array();
		$y = preg_match("/(".$beginstr.")(.*?)(".$endstr.")/s", $pagestring, $x);

		if(count($x)==4) {
			$selectedstr = $x[2];
			// get value
			$beginstr = "vulnerability \(";
			$endstr = "of 100";
			$val = array();
			$y = preg_match("/(".$beginstr.")(.*?)(".$endstr.")/s", $selectedstr, $val);			
			if(count($val) == 4) $myReturn['vulnerab'] = trim($val[2]);
			unset($val);
			
			// get reference
			$beginstr = "'>";
			$endstr = "</a>";
			$ref = array();
			$y = preg_match("#(".$beginstr.")(.*?)(".$endstr.")#i", $selectedstr, $ref);			
			if(count($ref) == 4) $myReturn['vulnerab_ref'] = trim($ref[2]);
			unset($ref);
		} 
		
		// TROPHIC LEVEL
		$beginstr = "<!-- Start resilience-->";
		$endstr = "<\/div>";
		$x = array();
		$y = preg_match("/(".$beginstr.")(.*?)(".$endstr.")/s", $pagestring, $x);

		if(count($x)==4) {
			$selectedstr = $x[2];
			// get value
			$beginstr = "doubling time ";
			$endstr = "\(";
			$val = array();
			$y = preg_match("/(".$beginstr.")(.*?)(".$endstr.")/s", $selectedstr, $val);			
			if(count($val) == 4) $myReturn['resil'] = trim($val[2]);
			unset($val);
			
			// get reference
			$beginstr = "\'>";
			$endstr = "<\/a>";
			$ref = array();
			$y = preg_match("#(".$beginstr.")(.*?)(".$endstr.")#i", $selectedstr, $ref);			
			if(count($ref) == 4) $myReturn['resil_ref'] = trim($ref[2]);
			unset($ref);
			
			// get K
			$beginstr = "K=";
			$endstr = ";";
			$k = array();
			$y = preg_match("/(".$beginstr.")(.*?)(".$endstr.")/s", $selectedstr, $k);			
			if(count($k) == 4) $myReturn['resil_K'] = trim($k[2]);
			
			// get tm
			$beginstr = "tm=";
			$endstr = ";";
			$tm = array();
			$y = preg_match("/(".$beginstr.")(.*?)(".$endstr.")/s", $selectedstr, $tm);			
			if(count($tm) == 4) $myReturn['resil_tm'] = trim($tm[2]);
			// get tmax
			$beginstr = "tmax=";
			$endstr = ";";
			$tmax = array();
			$y = preg_match("/(".$beginstr.")(.*?)(".$endstr.")/s", $selectedstr, $tmax);			
			if(count($tmax) == 4) $myReturn['resil_tmax'] = trim($tmax[2]);
			// get Fec
			$beginstr = "Fec=";
			$endstr = "\)";
			$fec = array();
			$y = preg_match("/(".$beginstr.")(.*?)(".$endstr.")/s", $selectedstr, $fec);			
			if(count($fec) == 4) $myReturn['resil_fec'] = trim($fec[2]);
		}
		
		// TROPHIC LEVEL
		$beginstr = "Trophic Level ";
		$endstr = "<\/div>";
		$x = array();
		$y = preg_match("/(".$beginstr.")(.*?)(".$endstr.")/s", $pagestring, $x);

		if(count($x)==4) {
			$selectedstr = $x[2];
			// get value
			$beginstr = "&nbsp;&nbsp;";
			$endstr = " &nbsp; ";
			$val = array();
			$y = preg_match("/(".$beginstr.")(.*?)(".$endstr.")/s", $selectedstr, $val);			
			if(count($val) == 4) $myReturn['tlevel'] = trim($val[2]);
			unset($val);
			
			// get standard error
			$beginstr = "&plusmn;";
			$endstr = "se;";
			$se = array();
			$y = preg_match("/(".$beginstr.")(.*?)(".$endstr.")/s", $selectedstr, $se);
			if(count($se) == 4) $myReturn['tlevelse'] = trim($se[2]);
			
			// get reference
			$beginstr = "\'>";
			$endstr = "<\/a>";
			$ref = array();
			$y = preg_match("#(".$beginstr.")(.*?)(".$endstr.")#i", $selectedstr, $ref);			
			if(count($ref) == 4) $myReturn['tlevelref'] = trim($ref[2]);
			unset($ref);
		}
		
		// PHYLOGENETIC DIVERSITY
		$beginstr = "Phylogenetic diversity index";
		$endstr = "<\/div>";
		$x = array();
		$y = preg_match("/(".$beginstr.")(.*?)(".$endstr.")/s", $pagestring, $x);

		if(count($x)==4) {
			$selectedstr = $x[2];
			//echo $selectedstr."\n";
			// get value
			$beginstr = "PD<sub>50<\/sub> = ";
			$endstr = "&nbsp;&nbsp;";
			$val = array();
			$y = preg_match("/(".$beginstr.")(.*?)(".$endstr.")/s", $selectedstr, $val);			
			if(count($val) == 4) $myReturn['pdiv'] = trim($val[2]);
			unset($val);
			
			// get reference
			$beginstr = "FBRefSummary.php\?ID=";
			$endstr = "'";
			$ref = array();
			$y = preg_match("/(".$beginstr.")(.*?)(".$endstr.")/s", $selectedstr, $ref);			
			if(count($ref) == 4) $myReturn['pdivref'] = trim($ref[2]);
			unset($ref);
		} 
		
		return($myReturn);
	}

	// extract lengths from string and return in units of mm & reference id & sex
	function getLen($sizeString) {
		$temp = explode(" ", $sizeString);
		
		$len = trim($temp[0]); 
		if(preg_match("/^[0-9,]+$/", $len)) $len = str_replace(",","",$len);
		
		
		if(trim($temp[1]) == "mm") $length = $len;
		else if(trim($temp[1]) == "cm") $length = $len * 10;
		else if(trim($temp[1]) == "m") $length = $len * 1000;
		else $length = 'units';
		
		$refId = getRef($sizeString);
		
		if(strpos($sizeString, "(male/unsexed,")!==false) $sex = 'male/unsexed';
		else if(strpos($sizeString, "(male,")!==false) $sex = 'male';
		else if(strpos($sizeString, "(female,")!==false) $sex = 'female';
		else $sex = NULL;
		
		return(array($length, $sex, $refId));
	}
	
	// extract weight from string and return in units of g & reference id
	function getWgt($wieghtString) {
		$wieghtString = str_replace("max. published weight:", "", $wieghtString);
		$temp = explode(" ", $wieghtString);
		
		$wgt = trim($temp[1]); 
		if(preg_match("/^[0-9,]+$/", $wgt)) $wgt = str_replace(",","",$wgt);
		
		if(trim($temp[2]) == "g") $weight = $wgt;
		else if(trim($temp[2]) == "kg") $weight = $wgt * 1000;
		else $weight = 'units';
		$refId = getRef($wieghtString);
		
		return(array($weight, $refId));
	}
?>