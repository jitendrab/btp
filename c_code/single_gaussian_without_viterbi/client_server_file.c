void writePlotFile(float **posterior, int totalNumFeatures, int *numStates){
  int                    i = 0, j = 0, s= 0, d = 0;
  int                    *featureState;
  FILE                   *plotFile;
  featureState = (int *)calloc(totalNumFeatures, sizeof(int *));
  for(i = 0; i < totalNumFeatures; i++){
    float max = -999999;
    int max_s = 0;
    for(s = 0; s < *numStates; s++){
      if(posterior[s][i] > max){
	max = posterior[s][i];
	max_s = s;
      }
    }
    featureState[i] = max_s;
  }
  plotFile = fopen("plot_data.txt", "w");
  if(!plotFile){
    fprintf(stderr, "unable to open plot file\n");
  }
  for(i = 0; i < totalNumFeatures; i++){
    fprintf(plotFile, "%d    %d\n", i, featureState[i]);
  }
  fclose(plotFile);
}

 writePlotFile(posterior, totalNumFeatures, numStates);
