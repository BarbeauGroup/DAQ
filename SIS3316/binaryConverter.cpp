//
//  binaryConverter.cpp
//  
//  Created by C. Awe on 5/14/2020
//  Based on code by grayson rich on 10/6/14.
//
//
//
//  modeled loosely after code by K. Kazkaz (LLNL)
//  original decoder by Kazkaz was used for 3320-generated data
//
//  this converter is meant to translate an NGM-binary file into a ROOT version
//  no analysis is performed, no data is removed except for NGM headers
//  The advantage here is that there is no need to compile, it can be run from a root session.
//
// USAGE: root[0] .L /path/to/this/file.cpp
//        root[1] convertData("/path/to/binary/file.bin")
//
#define DEBUG 0


//#include <stdio.h>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <stdio.h> 
#include <string.h> 

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TObject.h"
#include "TTreeIndex.h"

using namespace std;

//Supply this with the full path to your binary file and the path to where you'd like to store your output files.
void convertData( TString inputFilename, TString outputPath ){

    //Run specific DAQ settings.
    UShort_t numCards=1;
    UShort_t channelsPerCard=16;

    // check to see if the input file ends in '.bin' extension
    // if it doesn't, exit the program - user probably entered something wrong
    TString fileName_Core;
    TObjArray *inputParse1 = inputFilename.Tokenize('.');
    TString inputEnding = ((TObjString *) inputParse1->At(1))->GetString();
    if( inputEnding.EqualTo("bin") ){
	TString inputStart = ((TObjString *) inputParse1->At(0) )->GetString();
    	TObjArray *inputParse2 = inputStart.Tokenize('/');
	fileName_Core = ((TObjString *) inputParse2->At(inputParse2->GetLast()) )->GetString(); //Returns the "core" of the name i.e. "SIS3316Raw..."
    }
    else {
	cout << "Input file not the right format (doesn't end in .bin), quitting." << endl;
	exit(-1);
    }
    // generate the name for the output file
    outputPath += fileName_Core;
    outputPath += ".root";
    cout <<"The converted file will be: " << outputPath << endl;
    TFile *outFile = new TFile( outputPath,"RECREATE" );
        
    //Make ifstream object to read in file
    ifstream inFile;
    
    //Create output ttree variables
    ULong64_t timestamp;
    UShort_t peakHighIndex;
    UShort_t peakHighValue;
    UShort_t channelID;
    UChar_t formatBits;
    UInt_t accumulatorSum[6];
    UChar_t informationBits;
    UInt_t mawMaximumValue;
    UInt_t mawValueAfterTrigger;
    UInt_t mawValueBeforeTrigger;
    UInt_t mawTestData;
    UInt_t startEnergyValue;
    UInt_t maxEnergyValue;
    Bool_t mawTestFlag;
    Bool_t pileupFlag;
    UInt_t nSamples;
    UShort_t waveform[65536];
    
    //Create buffer object for reading in data
    UInt_t readBuffer[4096];
    char* bufferPointer = (char*)readBuffer;

    //Create variable to read data into from buffer
    UInt_t tmpWord;
    
    //Start by making a TTree of entries in order they arrive.
    //Then we'll have to sort them so they end up in chronological
    //order rather than channel/dump order
    TTree* unsortedTree = new TTree( "unsortedTree", "unsorted tree" );
    unsortedTree->SetDirectory(0);
    
    //Set branches
    unsortedTree->Branch( "channelID", &channelID, "channelID/s" );
    unsortedTree->Branch( "timestamp", &timestamp, "timestamp/l" );
    unsortedTree->Branch( "peakHighIndex", &peakHighIndex, "peakHighIndex/s" );
    unsortedTree->Branch( "peakHighValue", &peakHighValue, "peakHighValue/s" );
    unsortedTree->Branch( "formatBits", &formatBits, "formatBits/b" );
    unsortedTree->Branch( "accumulatorSum", accumulatorSum, "accumulatorSum[8]/i" );
    unsortedTree->Branch( "informationBits", &informationBits, "informationBits/b" );
    unsortedTree->Branch( "mawMaximumValue", &mawMaximumValue, "mawMaximumValue/i" );
    unsortedTree->Branch( "mawValueAfterTrigger", &mawValueAfterTrigger, "mawValueAfterTrigger/i" );
    unsortedTree->Branch( "mawValueBeforeTrigger", &mawValueBeforeTrigger, "mawValueBeforeTrigger/i" );
    unsortedTree->Branch( "mawTestData", &mawTestData, "mawTestData/i" );
    unsortedTree->Branch( "startEnergyValue", &startEnergyValue, "startEnergyValue/i" );
    unsortedTree->Branch( "maxEnergyValue", &maxEnergyValue, "maxEnergyValue/i" );
    unsortedTree->Branch( "mawTestFlag", &mawTestFlag, "mawTestFlag/O" );
    unsortedTree->Branch( "pileupFlag", &pileupFlag, "pileupFlag/O" );
    unsortedTree->Branch( "nSamples", &nSamples, "nSamples/i" );
    unsortedTree->Branch( "waveform", waveform, "waveform[nSamples]/s");

    //Actual tree we'll be writing, after sorting
    TTree* sis3316tree = (TTree*)unsortedTree->CloneTree( 0 );
    sis3316tree->SetName( "sis3316tree" );
    sis3316tree->SetDirectory( outFile );
    
    //	keep track of the time taken to process things
    // start a timer now
    Long64_t index = 0;
    time_t startTime, endTime;
    time_t processStartingTime;
    time( &processStartingTime );
    
    //packetWords stores number of expected word in dump
    UInt_t packetWords;
    //Determines how many event words we expect based on format bits
    UInt_t eventWordsFromFormatBits;
    
    //If there are less words than expected, set this flag to ignore
    //last bit of data. We have no idea of how ending data would be 
    //effected, depends on where in dump the crash happens
    Bool_t incompleteDumpFlag=0;

    //Dump number
    Long64_t spillNumber = 0;
    
    //Open file
    inFile.open( inputFilename.Data(), ifstream::in ); 


    //Get size of file, from:
    //https://stackoverflow.com/questions/2409504/using-c-filestreams-fstream-how-can-you-determine-the-size-of-a-file
    inFile.ignore(std::numeric_limits<std::streamsize>::max());
    std::streamsize bitsInFile = inFile.gcount();
    inFile.clear();
    inFile.seekg(0,std::ios_base::beg);

    printf( "Processing %s...\n", inputFilename.Data() );

    // seek past the header at the start of the binary file
    // 100 words = 400 bytes
    inFile.seekg( 400 );
    bitsInFile-=400;
    
    printf( "Header seeked through..\n" );
    
    while( inFile.good() ) {

      if (DEBUG) {
        std::cout<<"bits left in file: "<<bitsInFile<<endl;
      }
      
      // we're at the start of a spill
      // read 10 word spill header
      inFile.read( bufferPointer, 40 );
      bitsInFile-=40;
        
      if( DEBUG ) {
        // if we're debugging, print out the first 10 words of the first spill
        for( Int_t i = 0; i < 10; i++ ) {
          memcpy( &tmpWord, &readBuffer[i], 4 );
          printf( "Word %i of %llu spill: \t %x\n", i, spillNumber, tmpWord );
        }
      }
        
      //check if first word we read was EOF, shouldn't this be done prior to ABAB part? probably doesn't matter
      memcpy( &tmpWord, &readBuffer[0], 4 );
      if( tmpWord == 0x0E0F0E0F ) {
        std::cout<<"EOF reached"<<std::endl;
        break;
      }
        
      for (Int_t cardNumber = 0; cardNumber < numCards; cardNumber++) {
        
        // skip the packet header for the spill, one generated for each 3316 card
        inFile.read( bufferPointer, 8 );
        bitsInFile-=8;

        // now that we're inside of the spill, we have to parse data for each channel
        for( Int_t channelNumber = 0; channelNumber < channelsPerCard; channelNumber++ ) {
            
          //read packet header for channel, 8 words?
          inFile.read( bufferPointer, 32 );
          bitsInFile-=32;
			
          //8th word is the number of packet words for this channel
          memcpy( &tmpWord, &readBuffer[7], 4 );
          packetWords = tmpWord;
          bitsInFile-=(packetWords*4);
				
			
          if( DEBUG ) {
            std::cout<<"Bits in file before digitizer dump: "<<bitsInFile+packetWords*4<<std::endl;
            printf( "Number of words in packet for channel %i in spill %llu: \t %u\n", channelNumber, spillNumber, packetWords );
            std::cout<<"Bits in file after digitizer dump: "<<bitsInFile<<std::endl;
          }
			
          //Check if more words are expected than are in file. Minimum of 10 words expected for this channel, = 40bits
          //Set flag so we don't sort the channels we already have data for, i.e. throw out entire spill 
          if (bitsInFile < 40) {
            incompleteDumpFlag=1;
            std::cout<<"Found incomplete dump, remaining data not added to TTree"<<std::endl;
            break;
          }
				
          while( packetWords > 0 ) {
				
					  if( index % 100000 == 0 ) {
						  time(&endTime);
						  printf( "Processed %lli events in %i seconds\n", index, (Int_t)difftime(endTime, processStartingTime ));
					  }
				
					  // first two words of an event are there no matter what the format bits are set to
					  inFile.read( bufferPointer, 8 );
					  packetWords -= 2;
					  memcpy( &tmpWord, &readBuffer[0], 4 );
				  
					  formatBits = (UChar_t)(tmpWord & 0xf);
				  
					  channelID = (UShort_t)((tmpWord & 0xfff0) >> 4);
				  
					  timestamp = (static_cast<ULong64_t>(tmpWord & 0xffff0000)) << 16;
				  
					  memcpy( &tmpWord, &readBuffer[1], 4 );
				  
					  timestamp = timestamp | tmpWord;
				
					  if( DEBUG ) {
						  // print out the first two words of the event
						  memcpy( &tmpWord, &readBuffer[0], 4 );
						  printf("First two words of event:\t %x\t", tmpWord );
						  memcpy( &tmpWord, &readBuffer[1], 4 );
						  printf( "%x\n", tmpWord );
					  }
				  
					  if( DEBUG ) {
						  // print out the determined format bits
						  printf( "Format bits: %i\n", formatBits );
					  }
				
				
					  // now based on format bits we can determine the number of words to read per event
					  eventWordsFromFormatBits = 0;
					  if( (formatBits & 0x1) != 0 ) {
					  
						  if( DEBUG ) {
							  printf( "Reading words for format bit 0\n");
						  }
						  eventWordsFromFormatBits += 7;
						  inFile.read( bufferPointer, 7 * 4 );
						  packetWords -= 7;
						  memcpy( &tmpWord, &readBuffer[0], 4);
						  peakHighValue = tmpWord & 0xffff;
						  peakHighIndex = (tmpWord & 0xffff0000) >> 16;
						  if( DEBUG ) {
							  printf("Peak index: %i peak value: %i\n", peakHighIndex, peakHighValue );
						  }
						  memcpy( &tmpWord, &readBuffer[1], 4 );
						  informationBits = ( tmpWord & 0xff000000 ) >> 24;
						  accumulatorSum[0] = ( tmpWord & 0xffffff );
						  memcpy( &accumulatorSum[1], &readBuffer[2], 4 );
						  memcpy( &accumulatorSum[2], &readBuffer[3], 4 );
						  memcpy( &accumulatorSum[3], &readBuffer[4], 4 );
						  memcpy( &accumulatorSum[4], &readBuffer[5], 4 );
						  memcpy( &accumulatorSum[5], &readBuffer[6], 4 );
					  }
					  else {
						  peakHighIndex = 0;
						  peakHighValue = 0;
						  informationBits = 0;
						  accumulatorSum[0] = 0;
						  accumulatorSum[1] = 0;
						  accumulatorSum[2] = 0;
						  accumulatorSum[3] = 0;
						  accumulatorSum[4] = 0;
						  accumulatorSum[5] = 0;
					  }
					  if( (formatBits & 0x2) != 0 ) {
						  if( DEBUG ) {
							  printf( "Reading words for format bit 1\n");
						  }
						  eventWordsFromFormatBits += 2;
						  inFile.read( bufferPointer, 2 * 4 );
						  packetWords -= 2;
						  memcpy( &accumulatorSum[6], &readBuffer[0], 4 );
						  memcpy( &accumulatorSum[7], &readBuffer[1], 4 );
					  }
					  else {
						  // populate accumulators with 0 if they're not defined
						  accumulatorSum[6] = 0;
						  accumulatorSum[7] = 0;
					  }
					  if( (formatBits & 0x4) != 0 ) {
						  if( DEBUG ) {
							  printf( "Reading words for format bit 2\n");
						  }
						  eventWordsFromFormatBits += 3;
						  inFile.read( bufferPointer, 3 * 4 );
						  packetWords -= 3;
						  memcpy( &mawMaximumValue, &readBuffer[0], 4 );
						  memcpy( &mawValueAfterTrigger, &readBuffer[1], 4 );
						  memcpy( &mawValueBeforeTrigger, &readBuffer[2], 4 );
					  }
					  else {
						  mawMaximumValue = 0;
						  mawValueAfterTrigger = 0;
						  mawValueBeforeTrigger = 0;
					  }
					  if( (formatBits & 0x8) != 0 ) {
						  if( DEBUG ) {
							  printf( "Reading words for format bit 3\n");
						  }
						  eventWordsFromFormatBits += 2;
						  inFile.read( bufferPointer, 2 * 4 );
						  packetWords -= 2;
						  memcpy( &startEnergyValue, &readBuffer[0], 4 );
						  memcpy( &maxEnergyValue, &readBuffer[1], 4 );
					  }
					  else {
						  startEnergyValue = 0;
						  maxEnergyValue = 0;
					  }
				  
					  // the next word will determine the number of sample words we read	
					  inFile.read( (char*)&tmpWord, 4 );
					  packetWords -= 1;
				
				  
				  
					  nSamples = 2 * (tmpWord & 0x3ffffff);
				  
					  if( DEBUG ) {
						  printf( "Determined there are %i sample words\n", nSamples );
					  }
				  
				  
					  pileupFlag = (tmpWord & 0x4000000 ) >> 26;
					  mawTestFlag = ( tmpWord & 0x8000000 ) >> 27;
				  
            for( Int_t i = 0; i < (nSamples / 2 ); i++ ) {
                inFile.read( (char*)&tmpWord, 4 );
                packetWords -= 1;
                waveform[i*2] = tmpWord & 0xffff;
                waveform[i*2 + 1] = (tmpWord & 0xffff0000) >> 16;
            }
				
					  //                    inFile.read( (char*)&mawTestData, 4 );
					  //                    packetWords -= 1;
					  mawTestData = 0;
				  
					  unsortedTree->Fill();
					  index++;
				  } //End while packet words > 0 loop
			  } //End channel loop
      } //End digitizer card loop
        
        
        
      if (incompleteDumpFlag==false) {
		    // now we do the time ordering of the spill
		    if( DEBUG ) {
          printf( "Creating TTree index for spill %lli\n", spillNumber );
          time(&startTime);
		    }
		    unsortedTree->BuildIndex( "timestamp" );
		    TTreeIndex* treeIndex = (TTreeIndex*)unsortedTree->GetTreeIndex();
		    if( DEBUG ) {
          time(&endTime);
          printf( "Index created for spill %lli in %i seconds\n", spillNumber, (Int_t)difftime(endTime, startTime) );
		    }
		    
		    if( DEBUG ) {
          printf( "Tree index contains %lli entries\n", treeIndex->GetN() );
		    }
		    
		    for( Int_t i = 0; i < treeIndex->GetN(); i++ ) {
		        
          unsortedTree->GetEntry( treeIndex->GetIndex()[i] );
          if( DEBUG ) {
            printf( "%i-th event recalled is index %lli\n", i, treeIndex->GetIndex()[i] );
          }
          sis3316tree->Fill();   
		    } 
	    treeIndex->Delete();
      } //End incomplete dump flag if statement
      
      //If we had an incomplete dump for some reason, 
      //and it's NOT at the end of the file (what is going on with the digitizer?)
      //there's a chance there's still good data after? Reset
	    spillNumber++;
	    unsortedTree->Reset();
    } //End while inFile.good flag
    
    inFile.close();

    
    time( &endTime );
	cout << "Recorded " << index << " events in "
    << difftime( endTime, processStartingTime ) << " seconds" << endl << endl;
    
    unsortedTree->Delete();
    sis3316tree->Write("sis3316tree", TObject::kOverwrite);
    outFile->Close();

}

//This can be supplied with a .txt file listing multiple binary files. 
void processMultipleFiles( TString inputFilename, TString pathToOutputFiles ){

  //Step through the .txt file, processing as we go.
  ifstream file( inputFilename );
  string line;
  TString lineTString;
  Int_t lineNum = 0;
  while( file.good() ){
	if (getline ( file, line, '\n' )){
		lineTString = line;
		cout << "Processing file " << line << endl;
		convertData( lineTString, pathToOutputFiles );
	}
  }

  cout << "Routine complete, exiting." << endl;
}





