#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkCenteredEuler3DTransform.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkImageRegionIterator.h"
#include "itkCastImageFilter.h"

#include "itkSigmoidImageFilter.h"

#include "itkbasics.h"

#include "itkRayCastInterpolateImageFunction.h"

#include <boost/filesystem/operations.hpp>
#include <boost/format.hpp>


using std::string;


void usage()
{
  std::cerr << "\n";
  std::cerr << "Usage: DRR <options> [inputdir]\n";
  std::cerr << "       calculates the Digitally Reconstructed Radiograph from  \n";
  std::cerr << "       a volume. \n\n";
  std::cerr << "   where <options> is one or more of the following:\n\n";
  std::cerr << "       <-h>                    Display (this) usage information\n";
  std::cerr << "       <-v>                    Verbose output [default: no]\n";
  std::cerr << "       <-res float float>      Pixel spacing of the output image [default: 1x1mm]  \n";
  std::cerr << "       <-size int int>         Dimension of the output image [default: 501x501]  \n";
  std::cerr << "       <-sid float>            Distance of ray source (focal point) [default: 400mm]\n";
  std::cerr << "       <-t float float float>  Translation parameter of the camera \n";
  std::cerr << "       <-rx float>             Rotation around x,y,z axis in degrees \n";
  std::cerr << "       <-ry float>\n";
  std::cerr << "       <-rz float>\n";
  std::cerr << "       <-normal float float>   The 2D projection normal position [default: 0x0mm]\n";
  std::cerr << "       <-cor float float float> The centre of rotation relative to centre of volume\n";
  std::cerr << "       <-threshold float>      Threshold [default: 0]\n";
  std::cerr << "       <-o file>               Output image filename\n\n";
  std::cerr << "                               by  thomas@hartkens.de\n";
  std::cerr << "                               and john.hipwell@kcl.ac.uk (CISG London)\n\n";
  exit(1);
}

int main( int argc, char *argv[] )
{
  char *input_name = NULL;
  char *output_name = NULL;

  bool ok;
  bool verbose = false;

  float rx = 0.;
  float ry = 0.;
  float rz = 0.;

  float tx = 0.;
  float ty = 0.;
  float tz = 0.;

  float cx = 0.;
  float cy = 0.;
  float cz = 0.;

  float sid = 400.;

  float sx = 1.;
  float sy = 1.;

  int dx = 501;
  int dy = 501;

  float o2Dx = 0;
  float o2Dy = 0;

  double threshold=0;


  // Parse command line parameters

  while (argc > 1)
    {
    ok = false;

    if ((ok == false) && (strcmp(argv[1], "-h") == 0))
      {
      argc--; argv++;
      ok = true;
      usage();
      }

    if ((ok == false) && (strcmp(argv[1], "-v") == 0))
      {
      argc--; argv++;
      ok = true;
      verbose = true;
      }

    if ((ok == false) && (strcmp(argv[1], "-rx") == 0))
      {
      argc--; argv++;
      ok = true;
      rx=atof(argv[1]);
      argc--; argv++;
      }

    if ((ok == false) && (strcmp(argv[1], "-ry") == 0))
      {
      argc--; argv++;
      ok = true;
      ry=atof(argv[1]);
      argc--; argv++;
      }

    if ((ok == false) && (strcmp(argv[1], "-rz") == 0))
      {
      argc--; argv++;
      ok = true;
      rz=atof(argv[1]);
      argc--; argv++;
      }

    if ((ok == false) && (strcmp(argv[1], "-threshold") == 0))
      {
      argc--; argv++;
      ok = true;
      threshold=atof(argv[1]);
      argc--; argv++;
      }

    if ((ok == false) && (strcmp(argv[1], "-t") == 0))
      {
      argc--; argv++;
      ok = true;
      tx=atof(argv[1]);
      argc--; argv++;
      ty=atof(argv[1]);
      argc--; argv++;
      tz=atof(argv[1]);
      argc--; argv++;
      }

    if ((ok == false) && (strcmp(argv[1], "-cor") == 0))
      {
      argc--; argv++;
      ok = true;
      cx=atof(argv[1]);
      argc--; argv++;
      cy=atof(argv[1]);
      argc--; argv++;
      cz=atof(argv[1]);
      argc--; argv++;
      }

    if ((ok == false) && (strcmp(argv[1], "-res") == 0))
      {
      argc--; argv++;
      ok = true;
      sx=atof(argv[1]);
      argc--; argv++;
      sy=atof(argv[1]);
      argc--; argv++;
      }

    if ((ok == false) && (strcmp(argv[1], "-size") == 0))
      {
      argc--; argv++;
      ok = true;
      dx=atoi(argv[1]);
      argc--; argv++;
      dy=atoi(argv[1]);
      argc--; argv++;
      }

    if ((ok == false) && (strcmp(argv[1], "-sid") == 0))
      {
      argc--; argv++;
      ok = true;
      sid=atof(argv[1]);
      argc--; argv++;
      }

    if ((ok == false) && (strcmp(argv[1], "-normal") == 0))
      {
      argc--; argv++;
      ok = true;
      o2Dx=atof(argv[1]);
      argc--; argv++;
      o2Dy=atof(argv[1]);
      argc--; argv++;
      }

    if ((ok == false) && (strcmp(argv[1], "-o") == 0))
      {
      argc--; argv++;
      ok = true;
      output_name = argv[1];
      argc--; argv++;
      }

    if (ok == false) 
      {

      if (input_name == NULL) 
        {
        input_name = argv[1];
        argc--;
        argv++;
        }
      
      else 
        {
        std::cerr << "ERROR: Can not parse argument " << argv[1] << std::endl;
        usage();
        }
      }
    } 

  if (verbose) 
    {
    if (input_name)  std::cout << "Input image: "  << input_name  << std::endl;
    if (output_name) std::cout << "Output image: " << output_name << std::endl;
    }


  const     unsigned int   Dimension = 3;
  typedef short InputPixelType;
  typedef int FilterPixelType;
  typedef itk::Image< InputPixelType,  Dimension >   InputImageType;
  typedef itk::Image< FilterPixelType, Dimension >   FilterImageType;

  InputImageType::Pointer image;


  if (input_name) 
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
      } else if (boost::filesystem::is_directory( input_name )) {
        itkBasic::SeriesReader sreader( input_name );
	sreader.readSeriesData( 2 );
	try 
	{
	    itkBasic::ReaderType::Pointer imageReader = itkBasic::ReaderType::New();
	    itkBasic::FileNamesContainer fc;
	    sreader.getSeriesFileNames(0, fc);
	    image = itkBasic::getDicomSerie( fc, imageReader, 1 ); 
	}
	catch( itk::ExceptionObject & err ) 
	  { 
	  std::cerr << "ERROR: ExceptionObject caught !" << std::endl; 
	  std::cerr << err << std::endl; 
	  return EXIT_FAILURE;
	  } 
      } else {
	std::cerr << "input is neither dir nor file - please go away!" << std::endl;
	usage();
      }
    }
  
  // Print out the details of the input volume

  if (verbose) 
    {
    unsigned int i;
    const InputImageType::SpacingType spacing = image->GetSpacing();  
    std::cout << std::endl << "Input ";
    
    InputImageType::RegionType region = image->GetBufferedRegion();
    region.Print(std::cout);
    
    std::cout << "  Resolution: [";
    for (i=0; i<Dimension; i++) 
      {
      std::cout << spacing[i];
      if (i < Dimension-1) std::cout << ", ";
      }
    std::cout << "]" << std::endl;
    
    const InputImageType::PointType origin = image->GetOrigin();
    std::cout << "  Origin: [";
    for (i=0; i<Dimension; i++) 
      {
      std::cout << origin[i];
      if (i < Dimension-1) std::cout << ", ";
      }
    std::cout << "]" << std::endl<< std::endl;
    }

  typedef itk::ResampleImageFilter<FilterImageType, FilterImageType > FilterType;

  FilterType::Pointer filter = FilterType::New();

  typedef itk::SigmoidImageFilter< InputImageType, FilterImageType > SigmoidCasterType;
  SigmoidCasterType::Pointer sigmoidcaster = SigmoidCasterType::New();
  sigmoidcaster->SetInput( image );
  sigmoidcaster->SetOutputMaximum( 4000 );
  sigmoidcaster->SetOutputMinimum( 1000 );

  filter->SetInput( sigmoidcaster->GetOutput() );
  filter->SetDefaultPixelValue( -1000 );

  typedef itk::CenteredEuler3DTransform< double >  TransformType;

  TransformType::Pointer transform = TransformType::New();

  transform->SetComputeZYX(true);

  TransformType::OutputVectorType translation;

  translation[0] = tx;
  translation[1] = ty;
  translation[2] = tz;

  // constant for converting degrees into radians
  const double dtr = ( vcl_atan(1.0) * 4.0 ) / 180.0;

  transform->SetTranslation( translation );
  transform->SetRotation( dtr*rx, dtr*ry, dtr*rz );

  InputImageType::PointType   imOrigin = image->GetOrigin();
  InputImageType::SpacingType imRes    = image->GetSpacing();

  typedef InputImageType::RegionType     InputImageRegionType;
  typedef InputImageRegionType::SizeType InputImageSizeType;

  InputImageRegionType imRegion = image->GetBufferedRegion();
  InputImageSizeType   imSize   = imRegion.GetSize();

  imOrigin[0] += imRes[0] * static_cast<double>( imSize[0] ) / 2.0; 
  imOrigin[1] += imRes[1] * static_cast<double>( imSize[1] ) / 2.0; 
  imOrigin[2] += imRes[2] * static_cast<double>( imSize[2] ) / 2.0; 

  TransformType::InputPointType center;
  center[0] = cx + imOrigin[0];
  center[1] = cy + imOrigin[1];
  center[2] = cz + imOrigin[2];

  transform->SetCenter(center);

  if (verbose) 
    {
    std::cout << "Image size: "
              << imSize[0] << ", " << imSize[1] << ", " << imSize[2] << std::endl
              << "   resolution: "
              << imRes[0] << ", " << imRes[1] << ", " << imRes[2] << std::endl
              << "   origin: "
              << imOrigin[0] << ", " << imOrigin[1] << ", " << imOrigin[2] << std::endl
              << "   center: "
              << center[0] << ", " << center[1] << ", " << center[2] << std::endl
              << "Transform: " << transform << std::endl;
    }


  typedef itk::RayCastInterpolateImageFunction<FilterImageType,double> InterpolatorType;

  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  interpolator->SetTransform(transform);

  interpolator->SetThreshold(threshold);
  InterpolatorType::InputPointType focalpoint;

  focalpoint[0]= imOrigin[0];
  focalpoint[1]= imOrigin[1];
  focalpoint[2]= imOrigin[2] - sid/2.;

  interpolator->SetFocalPoint(focalpoint);

  if (verbose)
    {
    std::cout << "Focal Point: " 
              << focalpoint[0] << ", " 
              << focalpoint[1] << ", " 
              << focalpoint[2] << std::endl;
    }

  filter->SetInterpolator( interpolator );
  filter->SetTransform( transform );

  // setup the scene
  InputImageType::SizeType   size;

  size[0] = dx;  // number of pixels along X of the 2D DRR image 
  size[1] = dy;  // number of pixels along Y of the 2D DRR image 
  size[2] = 1;   // only one slice

  filter->SetSize( size );

  InputImageType::SpacingType spacing;

  spacing[0] = sx;  // pixel spacing along X of the 2D DRR image [mm]
  spacing[1] = sy;  // pixel spacing along Y of the 2D DRR image [mm]
  spacing[2] = 1.0; // slice thickness of the 2D DRR image [mm]

  filter->SetOutputSpacing( spacing );

  if (verbose) 
    {
    std::cout << "Output image size: " 
              << size[0] << ", " 
              << size[1] << ", " 
              << size[2] << std::endl;
    
    std::cout << "Output image spacing: " 
              << spacing[0] << ", " 
              << spacing[1] << ", " 
              << spacing[2] << std::endl;
    }

  double origin[ Dimension ];

  origin[0] = imOrigin[0] + o2Dx - sx*((double) dx)/2.; 
  origin[1] = imOrigin[1] + o2Dy - sy*((double) dy)/2.; 
  origin[2] = imOrigin[2] + sid/2.; 

  filter->SetOutputOrigin( origin );

  if (verbose) 
    {
    std::cout << "Output image origin: " 
              << origin[0] << ", " 
              << origin[1] << ", " 
              << origin[2] << std::endl;
    }

  // create writer

  if (output_name) 
    {
    FilterImageType::Pointer result = filter->GetOutput();
    
    sigmoidcaster->SetBeta( -500 );
    sigmoidcaster->SetAlpha( 5 );
    result->Update();
    itkBasic::ImageSave( itkBasic::ImageSharp(result, 0.1), boost::str( boost::format("%s.%s.png") % output_name % "lung" ), 1, 0); // Auto Level

    sigmoidcaster->SetBeta( 1000 );
    sigmoidcaster->SetAlpha( 300 );
    result->Update();
    itkBasic::ImageSave( itkBasic::ImageSharp(result, 0.1), boost::str( boost::format("%s.%s.png") % output_name % "bone" ), 1, 0); // Auto Level
    
    sigmoidcaster->SetBeta( 0 );
    sigmoidcaster->SetAlpha( 2000 );
    result->Update();
    itkBasic::ImageSave( itkBasic::ImageSharp(result, 0.1), boost::str( boost::format("%s.%s.png") % output_name % "normal" ), 1, 0); // Auto Level
    }
}
