struct stat sb;
bool pod_exist;
if (stat("./ITHACAoutput/POD", &sb) == 0 && S_ISDIR(sb.st_mode))
{
  pod_exist = true;
}else
{
        pod_exist = false;
        Info << "POD don't exist, performing a POD decomposition" << endl;
        mkDir("./ITHACAoutput/POD");

	// Creating a symbolic link for the folders:
	
        system("ln -s ../../constant ./ITHACAoutput/POD/constant");
        system("ln -s ../../0 ./ITHACAoutput/POD/0");
        system("ln -s ../../system ./ITHACAoutput/POD/system");
}
if(pod_exist == 1)
{
    Info << "The POD has already been performed, so just do the projection " <<endl;
}
