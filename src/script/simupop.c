// Fonctions

// Gère la sortie
function (void)afficheSortie(numeric alpha, numeric vecM, numeric vecE, numeric vecMG, numeric vecEG, numeric vecMF, numeric vecEF, numeric vecHistEnf)
{
	
	g = integerDiv(sim.generation, ENR)+1;
	
	cat(sim.generation + " " + alpha + " " + cor(vecMF,vecEF) +
	 " " + cor(vecMG,vecEG) + " " + cor(vecM,vecE) + " ");
	
	for (i in seqLen(size(vecHistEnf)))
	{
		cat(asString(vecHistEnf[i]) + " ");
	}
	catn();
	
}

// Chaque individu a une propriété "tagF" qui contient la taille de sa fratrie
// sauf à la première génération, où ces tagFs sont tirés au hasard dans une loi de Poisson

initialize() {
	
	initializeSLiMModelType("nonWF"); // non WF simulation
	initializeSLiMOptions(keepPedigrees=T); // keep pedigree to prevent inbreeding		
	initializeSex("A"); // autosome
	defineConstant("K", tailleBash);	// K = taille de pop
	defineConstant("ENR", enr_Bash); // nombre d'enregistrements
	defineConstant("TRANS", trans_Bash); // alpha
	defineConstant("PARENT", trans_parent); // parent who transmits
	defineConstant("DISTRIB", distrib_Bash); // distribution poisson ou gamma

	// autres constantes
	defineConstant("LAMBDA", 2); // paramètre pour loi de poisson

	// neutral mutations
	initializeMutationType("m1", 0.5, "f", 0.0);
	m1.convertToSubstitution = T;
	initializeTreeSeq(); // avec enregistrement d'arbre
	
	initializeGenomicElementType("g1", m1, 1.0);
	initializeGenomicElement(g1, 0, L_Bash); // taille du génome
	initializeMutationRate(0); // Pas de mutations, on les ajoute après sur l'arbre si besoin (plus rapide)
	initializeRecombinationRate(rec_Bash); // recombinaison
}

1 early() {
	// population
	sim.addSubpop("p1", K);
	// tagF : taille de fratrie - A la génération 1 : fixée à 1 pour tous (tester d'autres cas)
	p1.individuals.tagF = rep(1, p1.individualCount);
	// x : nombre d'enfants fait par l'individu (initialisé à 0)
	p1.individuals.x = rep(0, p1.individualCount);
}

reproduction(NULL, NULL) 
{
	
	self.active = 0;  // un seul appel de la fonction "reproduction" par génération


	//######## Section 1 : creation des couples ########//
	
	
	
	inds = p1.individuals;
    indsF = inds[inds.sex=='F'];
    indsM = inds[inds.sex=='M'];
    
    tailleF = length(indsF); // nb de femmes
	tailleM = length(indsM); // nb d'hommes
	
	
    // Calcul du alpha
    
    debut = deb_Bash; // mettre en constante
	fin = fin_Bash;
    
    if (sim.generation > debut & sim.generation <= fin)
	{
		alpha = trans_Bash;
	}
	else
	{
		alpha = 0;
	}	


	// pour avoir le nombre max de couples (inutile si sex-ratio = 0.5)
	
	if (tailleF < tailleM) 
	{
		taille = tailleF;
	}
	else
	{
		taille = tailleM;
	}
	
	
	indsF = p1.sampleIndividuals(taille, sex='F'); // On mélange (VERIFIER que ca melange)
	indsH = p1.sampleIndividuals(taille, sex='M');
	
	ex = p1.sampleIndividuals(1);
	couplesPop = matrix(rep(ex,taille*2), nrow=taille, ncol=2); // matrice pour les couples
	
	for (i in seqLen(taille))  // Formation des couples
	{
		couplesPop[i,] = c(indsF[i], indsM[i]);		
	}
	
	
	//### Section  2 : normalisation entre 0 et 1 pour calculer proba de se reproduire ###//
	
	proba = c();
	
	for (i in seqLen(taille))  // calcul proba de chaque couple non normalisee
	{		
		pere = couplesPop[i,1];
		mere = couplesPop[i,0];
		
		if (DISTRIB == 1) // if distribution is Poisson like
		{
			distrib = 1;
		}
		else
		{
			distrib = rgamma(1, 1, 1);
		}
		
		
		
		// Type de transmission
		
		if (PARENT == 0) // Transmission par la mère
		{
			p = distrib * ((mere.tagF)^(alpha));
		}
		else if (PARENT == 1) // Transmission par le père
		{
			p = distrib * ((pere.tagF)^(alpha));
		}
		else // Transmission biparentale
		{
			p = distrib * (((pere.tagF + mere.tagF)/2)^(alpha));
		}
		
		
		
		proba = c(proba, p); 		
	}
	
	// div vectorielle possible
	somme = sum(proba); // english
	
	probaNorm = proba / somme; // normalisation

	
	probaSum = cumSum(probaNorm); // somme cumulée des probas

	
	
	nb_children_done = 0;
	listeEnfants = c();
	
	//### Section 3 : création de K enfants ###//
	
	while (nb_children_done < K) // Création de K enfants
	{
		if (integerMod(nb_children_done, 2) == 0)
		{
			sexe = "M";
		}
		else
		{
			sexe = "F";
		}
				
		//### Section 3a : tirage aleatoire selon les probas calculees avant ###//
		
		u = runif(1,0,1);
		
		indice = 0;
		
		for (i in seqLen(taille))
		{
			if (probaSum[i] > u)
			{
				indice = i;
				break;
			}
			
		}
		
		//###### Section 3b : creation de l'enfant ######//
		
		
		pere = couplesPop[indice,1]; // on prend le pere
		mere = couplesPop[indice,0]; // on prend la mere
		
		enfant = subpop.addCrossed(mere, pere, sex=sexe); // creation de l'enfant
		enfant.x = 0.0; // l'enfant a 0 enfant au départ
		listeEnfants = c(listeEnfants, enfant); // on l'ajoute à la liste des enfants
		
		pere.x = pere.x + 1.0; // incremente x (nb d'enfants faits)
		mere.x = mere.x + 1.0;
		
		pere.tag = pere.pedigreeID;
		mere.tag = mere.pedigreeID;
		
		nb_children_done = nb_children_done + 1; // increment du compte d'enfants faits au total
		
	}
	
	

	for (i in seqLen(nb_children_done)) // Transmission aux enfants du succès repro des parents
	{
		idPere = listeEnfants[i].pedigreeParentIDs[1]; // monogamie donc on peut prendre le pere
		pere = p1.sampleIndividuals(1, tag=idPere); // on pioche le pere (sex = "M")
		listeEnfants[i].tagF = pere.x; // on transmet son nombre d'enfants faits	
	}

	
	

	// Enregistrement dans le tableau si cette génération doit etre enregistrée
	
	if((sim.generation > 2 & integerMod(sim.generation, ENR) == 0) | (sim.generation == 2)) 
	{
		
		// vecteurs totaux
		vecM = c();
		vecF = c();
		vecE = c();
		for (i in seqLen(size(inds)))
		{
			grPereAFait = inds[i].tagF;
			vecM = c(vecM, grPereAFait);
			vecE = c(vecE, inds[i].x);		
		}
		
		// vecteurs pere garçon
		
		vecMG = c();
		vecFG = c();
		vecEG = c();
		
		gars = indsM;
		
		for (i in seqLen(size(gars)))
		{
			grPereAFait = gars[i].tagF;
			vecMG = c(vecMG, grPereAFait);
			vecEG = c(vecEG, gars[i].x);	
		}
		
		
		// vecteurs pere fille
		
		vecMF = c();
		vecFF = c();
		vecEF = c();
		
		filles = indsF;
		
		for (i in seqLen(size(filles)))
		{
			grPereAFait = filles[i].tagF;
			vecMF = c(vecMF, grPereAFait);
			vecEF = c(vecEF, filles[i].x);			
		}
		
		// Compter les enfants par couple
	
		// vecteur de taille N
		nb = nb_Bash; // taille du vecteur qui contiendra la distribution du nombre d'enfants faits
		vecHistEnf = rep(0,nb);
		
		for (i in seqLen(taille))
		{
			pere = couplesPop[i,1];
			val = asInteger(pere.x[0]);
			if (val < nb)
			{
				vecHistEnf[val] = vecHistEnf[val] + 1;
			}
		}	
		
		afficheSortie(alpha, vecM, vecE, vecMG, vecEG, vecMF, vecEF, vecHistEnf); // on envoie en sortie les infos sur cette génération (pour traitement sous bash puis R)
	}
    
	
}



early() {
	inds = p1.individuals;
	adults = inds[inds.age > 0];
	juvs = inds[inds.age == 0];
	
	// non-overlapping generations : on tue les adultes
	adults.fitnessScaling = 0.0;
	
	
}


