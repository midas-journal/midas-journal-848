#include "itkMultiLocalCreasenessImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkInvertIntensityImageFilter.h"
#include "itkImage.h"

int main(int argc, char* argv [])
{
  
  if ( argc != 5 )
    {
    std::cerr << "Incorrect parameters: \nUsage:" << 
    argv[0] << "Input-image sigma rho Output-image\n"<< std::endl;
    return -1;
    }

  try 
    {
    typedef itk::Image<float,2>              float2dImg;
    typedef itk::ImageFileReader<float2dImg> ReaderFilter;
    ReaderFilter::Pointer reader = ReaderFilter::New();
    reader->SetFileName(argv[1]);
    reader->Update();
    typedef itk::MultiLocalCreasenessImageFilter<float2dImg,float2dImg>
      RidgeFilter;
    RidgeFilter::Pointer creases = RidgeFilter::New();
    creases->SetSigma(atof(argv[2]));
    creases->SetRho(atof(argv[3]));
    creases->SetInput(reader->GetOutput());
    creases->Update();
    typedef itk::ImageFileWriter<float2dImg> WriterFilter;
    WriterFilter::Pointer writer = WriterFilter::New();
    writer->SetFileName(argv[4]);
    writer->SetInput(creases->GetOutput());
    writer->Update();
    }
  catch( itk::ExceptionObject& e ) 
    {
    std::cerr << e.what() << std::endl;
    return -1;
    }
  catch( ... ) 
    {
    std::cerr << "Problema\n";
    return -2;
    }
  return 0;
}
