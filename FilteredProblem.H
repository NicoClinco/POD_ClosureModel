class FilteredProblem
{
	public:
	
	FilteredProblem();
	FilteredProblem
	(
         const PtrList<volVectorField>& _uorgfield,
         const PtrList<volScalarField>& _heorgfield,
	 const PtrList<volVectorField>& _umodes,
	 const PtrList<volScalarField>& _hemodes
	);
	FilteredProblem(int argc,char* argv[]);

	~FilteredProblem() {};

	// FOM 
	PtrList<volVectorField> UorgField; // Original
	PtrList<volVectorField> UfiField; // Filtered velocity

	PtrList<volScalarField> HeOrgField; // Original 
	PtrList<volScalarField> HefiField; // Filtered energy

	// POD	
	PtrList<volVectorField> Umodes;
	PtrList<volScalarField> Hemodes;

	int NUMODES;
	int NHEMODES;

	PtrList<volScalarField> nu_list; // Viscosity-fields

	scalar closureViscosity = 0.0;
	
	void FilteringStepU(); // Filter the original velocity);
	
	void FilteringStepHe(); // Filter the original Energy 
	
// 	void InitializeFilterViscosity(const scalar& visc);

	void WriteFilteredSolutions();
};

/* To review this constructor */ 

FilteredProblem::FilteredProblem():
	NUMODES(5),
	NHEMODES(5)
{
  Info << "Default nOfModes for velocity : " << NUMODES 
	  << endl;
  Info << "Default nOfModes for entalphy : " << NHEMODES
	  <<endl;
};


FilteredProblem::FilteredProblem
( 
   const PtrList<volVectorField>& _uorgfield,
   const PtrList<volScalarField>& _heorgfield,
   const PtrList<volVectorField>& _umodes,
   const PtrList<volScalarField>& _hemodes
):
  UorgField(_uorgfield),
  HeOrgField(_heorgfield),
  Umodes(_umodes),
  Hemodes(_hemodes)
{
 if(!Umodes.empty())
 {
   NUMODES = Umodes.size();
 }
 if(!Hemodes.empty())
 {
   NHEMODES = Hemodes.size();
 }
};

void
FilteredProblem::FilteringStepU()
{
  
   UfiField.setSize(UorgField.size()); 
   Eigen::MatrixXd UmassMatrix = ITHACAutilities::getMassMatrix( Umodes, NUMODES,true ); // <Phi_i,Phi_j>
   Eigen::MatrixXd SystemMatrix(NUMODES,NUMODES);
   Eigen::MatrixXd rhs(NUMODES,1);
   

   label snapI = 0;
   for( snapI = 0; snapI < UorgField.size(); snapI ++)
   {
	for ( label j = 0 ; j <NUMODES ; j++ )
	 {
	   for (label k = 0; k<NUMODES ; k++ )
	    {
		 
		 SystemMatrix(j,k) = UmassMatrix(j,k) 
	         -fvc::domainIntegrate( fvc::laplacian(closureViscosity*Umodes[k]) & Umodes[j]).value();
	    } // end loop Cols
	   rhs(j,0) = fvc::domainIntegrate( UorgField[snapI] & Umodes[j] ).value();
	 }// end loop Rows

       // Solve the problem with Eigen & obtain Temporal Coefficients
       
       
       Eigen::MatrixXd Coeffs = SystemMatrix.colPivHouseholderQr().solve(rhs);
        
        
        UfiField.set
	 (
	   snapI,
	   UorgField[snapI].clone()
	 );
	
	// Set to zero the internalfield:
	//
	UfiField[snapI].ref()*=dimensionedScalar("zero",dimless,0.0);

	
 		
       // Reconstruct-phase
       for( label coefInd = 0 ; coefInd < Coeffs.rows() ; coefInd++)
	  {
	    UfiField[snapI] += Umodes[coefInd] * Coeffs(coefInd,0); // Reconstruct the fields
       	  }// End local reconstruction
       
    }// end snapshot cycle
   Info << "End----filtering velocity----"<<" nTotalSnapShot "<< snapI<<endl;

}// end method


void
FilteredProblem::FilteringStepHe()
{
  // Creating the mass-matrix
   HefiField.setSize(HeOrgField.size());
   Eigen::MatrixXd HemassMatrix = ITHACAutilities::getMassMatrix( Hemodes, NHEMODES,true); // <Phi_i,Phi_j>
   Eigen::MatrixXd SystemMatrix(NHEMODES,NHEMODES);
   Eigen::MatrixXd rhs(NHEMODES,1);
   
   Info << " Filtering the Scalar fields "<<endl;
   label snapI = 0;
  for(snapI = 0; snapI < HeOrgField.size(); snapI ++)
   {
        for ( label j = 0 ; j <NHEMODES ; j++ )
         {
           for (label k = 0; k<NHEMODES ; k++ )
            {
		 
                 SystemMatrix(j,k) = HemassMatrix(j,k)
                 -fvc::domainIntegrate( fvc::laplacian(closureViscosity*Hemodes[k]) * Hemodes[j]).value();
            } // end loop Cols
           rhs(j,0) = fvc::domainIntegrate( HeOrgField[snapI] * Hemodes[j] ).value();
         }// end loop Rows

       // Solve the problem with Eigen & obtain Temporal Coefficients
       
       Eigen::MatrixXd Coeffs = SystemMatrix.colPivHouseholderQr().solve(rhs);

        HefiField.set
         (
           snapI,
           HeOrgField[snapI].clone()
         );
	// The boundary is not touched ( again )
        HefiField[snapI].ref()*=dimensionedScalar("zero",dimless,0.0);

	
       // Reconstruct-phase
       for( label coefInd = 0 ; coefInd < Coeffs.rows() ; coefInd++)
          {
            HefiField[snapI] += Hemodes[coefInd] * Coeffs(coefInd,0); // Reconstru
          }// End local reconstruction
     
    }// end snapshot cycle
 Info << "----End the Scalar fields filtering step---- "<<"Total HeSnap: "<< snapI <<endl;
}// end method 

/*
void FilteredProblem::InitializeFilterViscosity(const scalar& visc)
{
 
  nu_list.setSize(UorgField.size()); // set the size for viscosity

  for (label nSnap = 0 ; nSnap < nu_list.size() ; nSnap ++ )
   {
    // IMPORTANT : SET THE TIME NAME ! 
     nu_list.set
     (
       nSnap,
       new volScalarField
       ( 
	IOobject
	 (
	   "FilteredViscosity",
	   UorgField[nSnap].time().timeName(),
	   UorgField[nSnap].mesh(),
	   IOobject::NO_READ,
	   IOobject::NO_WRITE
	 ),     
      UorgField[nSnap].mesh(),
      dimensionedScalar("nu",dimViscosity,visc)
     )
    );

   };
 
}
*/

void FilteredProblem::WriteFilteredSolutions()
{
  mkDir("./ITHACAOUTPUT/filteredU");
  mkDir("./ITHACAOUTPUT/filteredHE");
  Info << "Writing filtered solution in :" << "ITHACAOUTPUT" << endl;

  for(label snapJ = 0 ; snapJ <  UfiField.size() ; snapJ ++ )
   {
	   Info << "Writing the filtered solution : " << snapJ << endl;
	   if ( !UfiField[snapJ].empty()){
	   ITHACAstream::exportSolution( UfiField[snapJ],UfiField[snapJ].time().times()[snapJ+2].name(),"./ITHACAOUTPUT/filteredU");
	   ITHACAstream::exportSolution( HefiField[snapJ], HefiField[snapJ].time().times()[snapJ+2].name(),"./ITHACAOUTPUT/filteredHE");
	   }
   }
  Info << "endl"<<endl;
  
}
