//
//  sample.h
//  mdaToolExample
//
//  Created by Daniel Beatty on 3/27/07.
//  Copyright 2007 Texas Tech University. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#ifdef __cplusplus
//	#include <deque>
	#include <dcgRenaissance/matrix.hpp>
	using namespace std;	
#endif


@interface sample : NSObject {
#ifdef __cplusplus
	deque<id> *sampleVector;
	deque<id> *mean;
	deque<id> *scatterMatrix;
#else
	void *sampleVector;
	void *mean;
	void *scatterMatrix;
#endif
	int numberOfSamples;	
}
/*
	Init With Samples: Conducts the following calculations 
	prepares the cached variables.
		Compute the sample mean
		Compute the number of samples
		Determine Normalized mean vector
		Compute Scatter Matrix for samples.
*/
-(sample) initWithSampleVector:(matrix) aSampleVector;
-(sample) initWithSampleVector:(matrix) aSampleVector
						overallMean:(vector)overallMean;

@end
