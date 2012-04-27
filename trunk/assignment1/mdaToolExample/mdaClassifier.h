//
//  mdaClassifier.h
//  mdaToolExample
//
//  Created by Daniel Beatty on 3/27/07.
//  Copyright 2007 Texas Tech University. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import "sample.h"

@interface mdaClassifier : NSObject {
	NSMutableArray *setOfClassesOfSamples;
	matrix *backgroundScatter;
	matrix *whitteningScatter;
	matrix *wDiscriminant;
	matrix *lamdbaDiscriminant;
}

/*
	
*/
-(mdaClassifier *) initWithClasses:(sample *);
-(mdaClassifier *) initWithMatrix: (matrix *) aMatrix
				classVector:(vector) aVector;

@end
