#include "itkResampleImageFilter.h"
#include "itkSigmoidImageFilter.h"

#include "itkbasics.h"

#include "itkAccumulateImageFilter.h"
#include <itkExtractImageFilter.h>
#include <itkCenteredRigid2DTransform.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <boost/filesystem/operations.hpp>
#include <boost/format.hpp>


using std::string;


int main( int argc, char *argv[] )
{
  string input_name;
  string output_dir;
  if (argc == 3) {
    input_name = argv[1];
    output_dir = argv[2];
  }

  const     unsigned int   Dimension = 3;
  const     unsigned int   OutDimension = 2;
  typedef short InputPixelType;
  typedef int FilterPixelType;
  typedef itk::Image< InputPixelType,  Dimension >   InputImageType;
  typedef itk::Image< FilterPixelType, Dimension >   FilterImageType;
  typedef itk::Image< FilterPixelType, OutDimension >   OutFilterImageType;

  InputImageType::Pointer image;
  itk::MetaDataDictionary dict;


  if (input_name.size() && output_dir.size()) 
    {
      if (boost::filesystem::is_regular_file( input_name )) {
	typedef itk::ImageFileReader< InputImageType >  ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( input_name );
	try 
	  { 
	  reader->Update();
	  } 
	catch( itk::ExceptionObject & err ) 
	  { 
	  std::cerr << "ERROR: ExceptionObject caught !" << std::endl; 
	  std::cerr << err << std::endl; 
	  return EXIT_FAILURE;
	  } 
	image = reader->GetOutput();
	dict = reader->GetMetaDataDictionary();
      } else if (boost::filesystem::is_directory( input_name )) {
        itkBasic::SeriesReader sreader( input_name );
	sreader.readSeriesData( 2 );
	try 
	{
	    itkBasic::ReaderType::Pointer imageReader = itkBasic::ReaderType::New();
	    itkBasic::FileNamesContainer fc;
	    sreader.getSeriesFileNames(0, fc);
	    image = itkBasic::getDicomSerie( fc, imageReader, 1 ); 
	    dict = *((*imageReader->GetMetaDataDictionaryArray())[0]);
	}
	catch( itk::ExceptionObject & err ) 
	  { 
	  std::cerr << "ERROR: ExceptionObject caught !" << std::endl; 
	  std::cerr << err << std::endl; 
	  return EXIT_FAILURE;
	  } 
      }
    }
    
    if (!image) {
	std::cerr << argv[0] << ": input output" << std::endl;
	exit(1);
    }
  

  typedef itk::SigmoidImageFilter< InputImageType, FilterImageType > SigmoidCasterType;
  SigmoidCasterType::Pointer sigmoidcaster = SigmoidCasterType::New();
  
  sigmoidcaster->SetInput( image );
  sigmoidcaster->SetOutputMaximum( 4000 );
  sigmoidcaster->SetOutputMinimum( 1000 );

  
  typedef itk::AccumulateImageFilter< FilterImageType, FilterImageType > AccumulateFilter;
  AccumulateFilter::Pointer accumulator = AccumulateFilter::New();
  accumulator->SetAccumulateDimension(1);
  accumulator->SetInput( sigmoidcaster->GetOutput() );

  typedef itk::ExtractImageFilter< FilterImageType, OutFilterImageType > ExtractFilter;
  ExtractFilter::Pointer extractor = ExtractFilter::New();
  extractor->SetInput( accumulator->GetOutput() );
  FilterImageType::Pointer accuOut = accumulator->GetOutput();
  accuOut->UpdateOutputInformation();
  FilterImageType::RegionType extractRegion = accuOut->GetLargestPossibleRegion();
  
  extractRegion.SetSize(1,0);
  
  extractor->SetExtractionRegion( extractRegion );

  typedef itk::ResampleImageFilter<OutFilterImageType, OutFilterImageType > ResampleFilter;
  ResampleFilter::Pointer resampler = ResampleFilter::New();
  resampler->SetInput( extractor->GetOutput() );
  
  typedef itk::BSplineInterpolateImageFunction< OutFilterImageType > InterpolatorType;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator->SetSplineOrder(3);
  
  resampler->SetInterpolator( interpolator );
  OutFilterImageType::Pointer exOut = extractor->GetOutput();
  exOut->UpdateOutputInformation();
  
  typedef itk::CenteredRigid2DTransform< double > TransformType;
  TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();
  OutFilterImageType::PointType exOutCenter = exOut->GetOrigin();
  exOutCenter[0] += (exOut->GetLargestPossibleRegion().GetSize()[0]-1) * exOut->GetSpacing()[0] *.5;
  exOutCenter[1] += (exOut->GetLargestPossibleRegion().GetSize()[1]-1) * exOut->GetSpacing()[1] *.5;
  transform->SetCenter( exOutCenter );
  transform->SetAngleInDegrees( 180 );
  resampler->SetTransform( transform );
  resampler->SetOutputParametersFromImage( exOut );

  OutFilterImageType::SpacingType resampleSpacing = exOut->GetSpacing();
  resampleSpacing.Fill( std::min( resampleSpacing[0], resampleSpacing[1] ) );
  OutFilterImageType::SizeType resampleSize;
  resampleSize[0] = exOut->GetLargestPossibleRegion().GetSize()[0] * exOut->GetSpacing()[0] / resampleSpacing[0];
  resampleSize[1] = exOut->GetLargestPossibleRegion().GetSize()[1] * exOut->GetSpacing()[1] / resampleSpacing[1];
  resampler->SetSize( resampleSize );
  resampler->SetOutputSpacing( resampleSpacing );
  
  OutFilterImageType::Pointer result = resampler->GetOutput();
  
  sigmoidcaster->SetBeta( -500 );
  sigmoidcaster->SetAlpha( 5 );
  result->Update();

  int outDicomIndex = 0;
  itk::EncapsulateMetaData( dict, "0008|0008", string("DERIVED\\SECONDARY\\AXIAL"));
  
  boost::filesystem::path outpath = output_dir;
  outpath = outpath / "IM%06d";
  
  std::vector< itk::MetaDataDictionary* > dictArray;
  dictArray.push_back(&dict);
  
  itkBasic::writeDicomSeries( itkBasic::ImageRescale(itkBasic::ImageSharp(result, 0.5), -1000, 4000), outpath.string(), &dictArray, outDicomIndex);
//  itkBasic::ImageSave( itkBasic::ImageSharp(result, 0.5), boost::str( boost::format("%s.%s.png") % output_name % "lung" ), 1, 0); // Auto Level

  sigmoidcaster->SetBeta( 1000 );
  sigmoidcaster->SetAlpha( 300 );
  result->Update();
  itkBasic::writeDicomSeries( itkBasic::ImageRescale(itkBasic::ImageSharp(result, 0.5), -1000, 4000), outpath.string(), &dictArray, outDicomIndex);
//  itkBasic::ImageSave( itkBasic::ImageSharp(result, 0.5), boost::str( boost::format("%s.%s.png") % output_name % "bone" ), 1, 0); // Auto Level
  
  sigmoidcaster->SetBeta( 0 );
  sigmoidcaster->SetAlpha( 2000 );
  result->Update();
  itkBasic::writeDicomSeries( itkBasic::ImageRescale(itkBasic::ImageSharp(result, 0.5), -1000, 4000), outpath.string(), &dictArray, outDicomIndex);
//  itkBasic::ImageSave( itkBasic::ImageSharp(result, 0.5), boost::str( boost::format("%s.%s.png") % output_name % "normal" ), 1, 0); // Auto Level
}
