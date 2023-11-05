<?php
	# Written by Noel A. Heim, Tufts University, noel.heim@tufts.edu
	# Last modified 20 November 2020
	# Writes data to tab-delimited file called fishbase_data.txt
	
	include('fbFunctions.php');
	date_default_timezone_set("America/New_York");	

	$t0 = microtime(true);
	echo "FishBase extraction started: ".date("Y-m-d g:i A",$t0)."\n";
	
	$append_existing_file = true;
	$fishbase_max_id = 150001;
	$fishbase_min_id = 2;
	
	if($append_existing_file == true) {
		// get last id & set to min search id
		$fh2 = fopen("fishbase_data.txt", "r");
		while ($file_row = fgetcsv($fh2, 1000, "\t")) $fbid = $file_row[0];
        fclose($fh2);
		$fishbase_min_id = $fbid + 1;
		$fh = fopen("fishbase_data.txt", "a");
	} else {
		// set up new file
		$column_names = array(
			'fb_taxon_id','taxon_name','taxon_uri','phylum','class','order','family','genus',
			'sl_mm','sl_sex','sl_ref','tl_mm','tl_sex','tl_ref','fl_mm','fl_sex','fl_ref','dw_mm','dw_sex','dw_ref','max_weight_g','max_weight_ref',
			'benthic','pelagic','benthopelagic','bathypelagic','pelagic-neritic','pelagic-oceanic','demersal','reef-associated',
			'anadromous','catadromous','amphidromous','potamodromous','limnodromous','oceanodromous','non-migratory',
			'freshwater','brackish','marine',
			'use_fisheries','use_aquaculture','use_gamefish','use_bait','use_aquarium','use_ref',
			'phylo_div','phylo_div_ref',
			'troph_level','troph_level_se','troph_level_ref',
			'resil','resil_ref','resil_K','resil_tm','resil_tmax','resil_fec',
			'vulnerab','vulnerab_ref','price','price_ref',
			'main_ref'
		);
		$fh = fopen("fishbase_data.txt", "w");
		$file_row = implode("\t",$column_names)."\n";
		fwrite($fh, $file_row);	
	}
	// A file to store missed taxa
	$fh2 = fopen("missedTaxa.txt", "w");
	fwrite($fh2, "fb_id\tspecies_name\n");	
	
	# short code block used in debugging
	/* 
	$bluefintuna = 147;
	$clownfish = 9209;
	$sardine = 1043;
	$cod = 69;
	
	$randsample = array($bluefintuna,$clownfish,$sardine,$cod,28450,21336,57723,33921,45974,46860,13127,48328,63425,50937,30204,40670,57985,56968,8469,9460,34271,33060,37051,14862,4047,8724,60824,26034,55175,50589,16439,4347,55202,33283,23385,9477,14150,59247,508,44000,31616,62468,56767,9152,59134,41846,32873,36165,1998,15559,35206,19437,19755,22923,27070,63209,41675,20517,53285,56413,56208,16093,50149,32800,16546,36707,53533,27413,47996,48224,4826,35425,12682,18693,29289,2723,45425,29177,51672,63892,6240,13231,28716,11322,7628,20356,7953,11484,31247,7870,55243,36130,29187,5259,43655,52439,48244,6392,40300,63487,17044,51752,27728,15802,53381,44118,44725,29801,43774,22275,25435,7092,45846,52840,53635,2932,5201,29593,36072,6273,62716,37304,54018,14874,21440,38998,3687,35517,28244,26843,29125,62220,59342,313,24318,5211,36173,13490,49089,7157,8571,20433,15635,38498,32245,38042,44535,32438,43634,32170,5876,31465,24264,36652,14810,55656,20266,15820,27968,32906,35480,49721,45471,6297,15695,48338,29909,34905,59510,41472,41130,61067,18279,4200,53785,34411,63180,3966,4897,48194,29469,54426,47526,24084,45475,62845,20896,20527,10802,28786,9757,1234,10037,35388,21787,31146,37803,53728,21447,62434,51174,28811,39590,15870,33904,9117,26538,49544,2003,39044,28854,33943,20824,20699,29049,47567,55228,9480,48414,25549,15500,2240,51227,47623,27677,48081,19977,15880,29605,25305,16828,648,38468,55095,18973,5784,23175,23788,62025,1313,6887,31464,54942,27788,56128,32321,56062,34917,11369,36185,59120,43197,24566,6618,24377,2574,28747,6058,1017,13605,28218,12794,17271,21058,1474,14320,29215,35056,19075,26474,20875,33069,46789,12690,9839,43726,14313,28099,4808,46367,50011,37435,53209,4901,13143,36951,18654,34245,10463,62879,21742,13964,62123,28683,14388,51327,29067,16391,47550,59765,3860,59960,16351,22348,6656,24241,24516,50935,6246,17970,33417,1557,40340,21209,33011,46025,20988,17712,583,15062,32443,61568,18817,52687,13542,52152,57361,6581,55349,22365,42946,61182,20395,43804,31874,2289,42257,58052,51409,18091,39197,26733,60657,24745,43467,28309,14083,54842,53578,29791,59296,49720,33897,10978,42857,8988,32161,26775,34984,14621,62397,21679,41815,29653,19120,5795,33438,56296,8777,53257,13349,10918,36005,21029,19898,21073,13335,47436,14499,23160,37206,57759,29803,40425,24355,49485,63926,11270,35541,24864,58250,53302,49366,18752,4290,46271,41889,10923,40600,45537,17072,20033,20401,8526,19581,35679,14360,47205,25071,21471,42424,51155,52215,4359,16362,18528,26568,7113,11312,4498,9391,29345,30099,35845,23231,22707,53176,51627,16330,11685,58986,61539,62746,17531,63144,62968,53663,42074,23988,24469,5236,31983,11098,16008,55769,61808,13953,62137,23653,30094,24,23750,61095,17832,25011,252,48945,30126,40781,36371,52923,31057,38502,50948,56398,35720,57258,47254,15438,8088,58028,62457,33292,44828,61400,43583,15599,56721,53228,1374,56312,17769,60728,19169,5410,30377,28375,35176,26439,6771,27212,47080,21444,28409,40985,27655,14969,42235,37779,32729,25205,14304,48374,41373,15910,6092,61402,16306,60692,22995,28518,14071,1591,35193,55065,19172,21398,60149,4777,4822,7083,45775,53129,3952,17076,58076,11258,16485,31506,21180,21549,58286,5385,42514,15743,25882,26189,23699,57363,31509,5662,24170,15761,52281,31355,50753,14482,7032,15025,27251,759,36380,39014,24811,13644,19729,8676,25135,6012,23141,58219,31527,27029,21892,7808,40013,53900,27603,63880,63213,11401,52477,25997,47230,25895,26605,47894,23651,15811,40341,46090,18857,57463,63090,9365,16841,52482,41398,38177,28763,55273,17034,112,37853,62945,56512,23228,41781,37187,39783,13176,21042,41723,46339,57745,61562,13278,42286,50569,58796,3932,17678,6421,49723,13135,54839,29432,62402,33062,9519,57916,23117,5279,23055,34963,15330,38860,40261,7222,9118,22841,1829,50830,44697,49478,6804,21975,21556,42075,35670,2217,30506,26338,53514,48691,21795,964,12507,20070,60988,40138,34721,19547,14202,11441,19096,50534,58191,50441,18772,49507,8703,31123,49533,10630,38419,59734,38637,33727,59694,254,15264,38007,50150,47258,23208,1948,35748,38362,63475,2071,21401,165,54778,56246,18610,9013,58407,24678,57713,9478,8177,21254,53117,1055,17894,34396,10273,22336,15983,43403,51559,335,39556,17624,40675,59408,62316,35523,37985,9334,30255,55232,37944,43891,8272,20522,34597,38262,49105,29396,7065,21247,14395,37427,12982,37387,34322,7710,44,53613,49317,30486,60591,36232,43871,19453,14263,17157,19949,51697,34609,51005,62901,53269,2348,33158,2495,61140,9957,19642,754,10982,14820,9553,21510,63895,36141,27605,56834,50269,9633,39825,1406,13448,34814,7758,22425,1201,13410,57935,28339,2326,36310,52863,42669,56796,3982,11992,53354,50813,3117,13806,62887,22215,35830,8058,20711,58764,5576,57598,31273,44148,30623,6073,17991,23532,22741,31879,63565,14567,5382,44022,34391,29512,33199,30960,58499,13124,49063,13377,13194,38687,43358,11262,28359,41089,60537,59885,1568,34292,13214,40280,39140,22260,51505,52841,13109,47942,22976,57931,15634,26715,36841,16285,56070,53072,53373,21004,1827,52713,48771,57218,46021,18605,54004,23923,23448,38798,47696,2674,12449,32680,35966,28561,19018,21887,2827,33443,40754,17930,38033,34572,20685,39791,49786,31731,42924,32907,46900,52172,4204,44490,47972,15780,42040,42847,5191,42447,27278,33281,29904,45192,56104,59111,29117,48008,44287,15613,24083,1646,9269,4059,3278,58652,14319,54599,18628,27117,60494,16793,44293,62923,15435,16201,12431,35389,15763,18208,58105,52963,8321,53338,10751,5377,13499,46746,63639,20155,24117,45875,42885,29707,46273,37040,12932,42826,48697,35681,36600,43829,11382,32835,59426,4065,54049,16443,17105,57706,5179,11420,36390,34288,21783,50068,18282,5162,11355,31389,43898,57163,3849,55377,24159,62207,12170,42494,36776,24356,40699,48989,28429,5960,15714,1052,49997,9668,20580,56311,2547,50662,52627,29267,54740,15191,32974,15304,34578,10389,31827,50026,25311,26193,45034,25016,59767,45531,41688,24871,44393,20068,38183,19919,7452,55689,43147,7779,30292,5551,32862,26889,33826,59395,2065);
	$randsample = array($bluefintuna,$clownfish,$sardine,$cod);
	asort($randsample); echo "There are ".count($randsample)." randomly selected species.\n";
	*/
	
	for($i = $fishbase_min_id; $i < $fishbase_max_id; $i++) {				
		if($i % 500 == 0) echo $i."\n";
		
		libxml_use_internal_errors(true);
		$feed = simplexml_load_file("http://www.fishbase.us/webservice/summary/fb/showXML.php?identifier=FB-".$i);
		
		if(is_countable($feed) && count($feed) > 0 && count($feed->taxon) > 0) {		
			$taxon_uri = NULL;
			$size = NULL;
			$size_ref = NULL;
			$habitat = NULL;
			
			// taxon info from Darwin Core
			$dwc = $feed->taxon->children("http://rs.tdwg.org/dwc/dwcore/");
					
			$name = $dwc->ScientificName;
			if(isset($dwc->Phylum)) $phylum = $dwc->Phylum;
			else $phylum = NULL;
	
			if(isset($dwc->Class)) $class = $dwc->Class;
			else $class = NULL;
	
			if(isset($dwc->Order)) $order = $dwc->Order;
			else $order = NULL;
	
			if(isset($dwc->Family)) $family = $dwc->Family;
			else $family = NULL;
	
			if(isset($dwc->Genus)) $genus = $dwc->Genus;
			else $genus = NULL;
			
			// gen source uri from dc
			$dc = $feed->taxon->children("https://purl.org/dc/elements/1.1/");
			$taxon_uri = $dc->source;
			
			// get sizes & habitats
			$tl = NULL;
			$tl_ref = NULL;
			$tl_sex = NULL;
			$sl = NULL;
			$sl_ref = NULL;
			$sl_sex = NULL;
			$fl = NULL;
			$fl_ref = NULL;
			$fl_sex = NULL;
			$dw = NULL;
			$dw_ref = NULL;
			$dw_sex = NULL;
			$weight = NULL;
			$weight_ref = NULL;
			$benthic = NULL;
			$pelagic = NULL;
			$benthopelagic = NULL;
			$bathypelagic = NULL;
			$pelagic_neritic = NULL;
			$pelagic_oceanic = NULL;
			$demersal = NULL;
			$reef_assoc = NULL;
			$anadromous = NULL;
			$catadromous = NULL;
			$amphidromous = NULL;
			$potamodromous = NULL;
			$limnodromous = NULL;
			$oceanodromous = NULL;
			$non_migratory = NULL;
			$freshwater = NULL;
			$brackish = NULL;
			$marine = NULL;
			$phylo_div = NULL;
			$phylo_div_ref = NULL;
			$fisheries = NULL;
			$aquaculture = NULL;
			$gamefish = NULL;
			$bait = NULL;
			$aquarium = NULL;
			$use_ref = NULL;
	
			foreach($feed->taxon->dataObject as $do) {
				// Size
				if($do->children("https://purl.org/dc/elements/1.1/")->identifier == "FB-Size-".$i) {
					$size = $do->children("https://purl.org/dc/elements/1.1/")->description;
					$size = str_replace("unsexed;", "unsexed,", $size);
					$size = str_replace("male;", "male,", $size);
					$size = explode("; ", $size);
					foreach($size as $sizepart) {
						if(strpos($sizepart, " TL ")!==false) {
							$temp = getLen($sizepart);
							$tl = $temp[0];
							$tl_sex = $temp[1];
							$tl_ref = $temp[2];							
						} else if(strpos($sizepart, " SL ")!==false) {
							$temp = getLen($sizepart);
							$sl = $temp[0];
							$sl_sex = $temp[1];
							$sl_ref = $temp[2];
						} else if(strpos($sizepart, " FL ")!==false) {
							$temp = getLen($sizepart);
							$fl = $temp[0];
							$fl_sex = $temp[1];
							$fl_ref = $temp[2];
						} else if(strpos($sizepart, " DW ")!==false) {
							$temp = getLen($sizepart);
							$dw = $temp[0];
							$dw_sex = $temp[1];
							$dw_ref = $temp[2];
						} else if(strpos($sizepart, "max. published weight")!==false) {
							$temp = getWgt($sizepart);
							$weight = $temp[0];
							$weight_ref = $temp[1];
						}
					}
				}
				// Habitat
				if($do->children("https://purl.org/dc/elements/1.1/")->identifier == "FB-Habitat-".$i) {
					$habitat = strtolower($do->children("https://purl.org/dc/elements/1.1/")->description);

					if(preg_match("/\bbenthopelagic\b/i", $habitat)==1) $benthopelagic = true;
					if(preg_match("/\bbathypelagic\b/i", $habitat)==1) $bathypelagic = true;
					if(preg_match("/\bpelagic-neritic\b/i", $habitat)==1) $pelagic_neritic = true;
					if(preg_match("/\bpelagic-oceanic\b/i", $habitat)==1) $pelagic_oceanic = true;
					if(preg_match("/\bpelagic\b/i", $habitat)==1) $pelagic = true;
					if(preg_match("/\bbenthic\b/i", $habitat)==1) $benthic = true;
					if(preg_match("/\bdemersal\b/i", $habitat)==1) $demersal = true;
					if(preg_match("/\breef-associated\b/i", $habitat)==1) $reef_assoc = true;
					
					if(preg_match("/\banadromous\b/i", $habitat)==1) $anadromous = true;
					if(preg_match("/\bcatadromous\b/i", $habitat)==1) $catadromous = true;
					if(preg_match("/\bamphidromous\b/i", $habitat)==1) $amphidromous = true;
					if(preg_match("/\bpotamodromous\b/i", $habitat)==1) $potamodromous = true;
					if(preg_match("/\blimnodromous\b/i", $habitat)==1) $limnodromous = true;
					if(preg_match("/\boceanodromous\b/i", $habitat)==1) $oceanodromous = true;
					if(preg_match("/\bnon-migratory\b/i", $habitat)==1) $non_migratory = true;
					
					if(preg_match("/\bfreshwater\b/i", $habitat)==1) $freshwater = true;
					if(preg_match("/\bbrackish\b/i", $habitat)==1) $brackish = true;
					if(preg_match("/\bmarine\b/i", $habitat)==1) $marine = true;
				}
				
				// Fisheries
				if($do->children("https://purl.org/dc/elements/1.1/")->identifier == "FB-Uses-".$i) {
					$use_str = strtolower($do->children("https://purl.org/dc/elements/1.1/")->description);
					$ind_uses = explode(";", $use_str);
					foreach($ind_uses as $use) {					
						if(preg_match("/fisheries:/i", $use)==1) $fisheries = trim(ltrim(trim($use), "fisheries:"));
						if(preg_match("/aquaculture:/i", $use)==1) $aquaculture = trim(ltrim(trim($use), "aquaculture:"));
						if(preg_match("/gamefish:/i", $use)==1) $gamefish = trim(ltrim(trim($use), "gamefish:"));
						if(preg_match("/bait:/i", $use)==1) $bait = trim(ltrim(trim($use), "bait:"));
						if(preg_match("/aquarium:/i", $use)==1) $aquarium = trim(ltrim(trim($use), "aquarium:"));
					}
					$use_ref_str = $do->reference;
					foreach($use_ref_str as $refPart) $use_ref .= getRef($refPart).";";
					$use_ref = rtrim($use_ref, ";");
				}
				
			}			
			$alt_uri = "https://www.fishbase.de/summary/speciessummary.php?id=".$i;
			$summary_page = file_get_contents($alt_uri);
			if($summary_page === FALSE) fwrite($fh2, $i."\t".$name."\n");

			// get phylogenetic diversity, trophic level, resiliance, and vulnerability 
			$otherData = getPhyloDiv($summary_page, $i);
			
			$taxon_info = $i."\t".$name."\t".$taxon_uri."\t".$phylum."\t".$class."\t".$order."\t".$family."\t".$genus;
			$size_info = $sl."\t".$sl_sex."\t".$sl_ref."\t".$tl."\t".$tl_sex."\t".$tl_ref."\t".$fl."\t".$fl_sex."\t".$fl_ref."\t".$dw."\t".$dw_sex."\t".$dw_ref."\t".$weight."\t".$weight_ref;
			$habitat_info1 = $benthic."\t".$pelagic."\t".$benthopelagic."\t".$bathypelagic."\t".$pelagic_neritic."\t".$pelagic_oceanic."\t".$demersal."\t".$reef_assoc;
			$habitat_info2 = $anadromous."\t".$catadromous."\t".$amphidromous."\t".$potamodromous."\t".$limnodromous."\t".$oceanodromous."\t".$non_migratory;
			$habitat_info3 = $freshwater."\t".$brackish."\t".$marine;
			$uses = $fisheries."\t".$aquaculture."\t".$gamefish."\t".$bait."\t".$aquarium."\t".$use_ref;
			$other_info = implode("\t",$otherData);
			$file_row = $taxon_info."\t".$size_info."\t".$habitat_info1."\t".$habitat_info2."\t".$habitat_info3."\t".$uses."\t".$other_info."\n";

			fwrite($fh, $file_row);
		}
	}
	fclose($fh);
	fclose($fh2);

	// timer		
	$t1 = microtime(true);
	echo "FishBase scrape ended: ".date("Y-m-d g:i A",$t1)."\n";
	echo "\nThis process took ".time_elapsed($t1, $t0).".\n";
	
	// the end
	exit();	
?>