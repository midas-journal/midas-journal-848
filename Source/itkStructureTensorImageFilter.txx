#ifndef __itkStructureTensorImageFilter_txx
#define __itkStructureTensorImageFilter_txx

#include "itkStructureTensorImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "itkImageFileWriter.h"

namespace itk
{

/**
 * Constructor
 */
template <typename TInputImage, typename TOutputImage, typename TOutputImage2 >
StructureTensorImageFilter<TInputImage,TOutputImage,TOutputImage2>
::StructureTensorImageFilter()
{
  m_NormalizeAcrossScale = false;

  if( ImageDimension > 1)
    {
    m_SmoothingFilters.resize(ImageDimension);
    }

  size_t order[ImageDimension];  
  for( unsigned int i = 0; i<ImageDimension; ++i )
    {
    for( unsigned int j = 0; j<ImageDimension; ++j )
      order[j] = (j == i) ? 1 : 0;
    m_SmoothingFilters[ i ] = GaussianFilterType::New();
    m_SmoothingFilters[ i ]->SetOrder(order);
    }
  
  this->SetSigma( 1.0 );
  this->SetRho( 1.0 );
  
  this->SetNthOutput( 0, this->MakeOutput( 0 ) );
  this->SetNthOutput( 1, this->MakeOutput( 1 ) );
}

/**
 * Set value of Sigma
 */
template <typename TInputImage, typename TOutputImage, typename TOutputImage2 >
void
StructureTensorImageFilter<TInputImage,TOutputImage,TOutputImage2>
::SetSigma( RealType sigma )
{
  m_sigma = sigma;

  for( unsigned int i = 0; i<ImageDimension; i++ )
    {
    m_SmoothingFilters[ i ]->SetVariance( m_sigma );
    }
  this->Modified();
}

/**
 * Set value of Rho
 */
template <typename TInputImage, typename TOutputImage, typename TOutputImage2 >
void
StructureTensorImageFilter<TInputImage,TOutputImage,TOutputImage2>
::SetRho( RealType rho )
{
  m_rho = rho;
  this->Modified();
}

/**
 * Set Normalize Across Scale Space
 */
template <typename TInputImage, typename TOutputImage, typename TOutputImage2 >
void
StructureTensorImageFilter<TInputImage,TOutputImage,TOutputImage2>
::SetNormalizeAcrossScale( bool normalize )
{
  m_NormalizeAcrossScale = normalize;

  for( unsigned int i = 0; i<ImageDimension; i++ )
    {
    m_Smoothin/gFilters[ i ]->SetNormalizeAcrossScale( normalize );
    }
  this->Modified();

}

template <typename TInputImage, typename TOutputImage, typename TOutputImage2 >
void
StructureTensorImageFilter<TInputImage,TOutputImage,TOutputImage2>
::GenerateInputRequestedRegion() throw(InvalidRequestedRegionError)
{
  // call the superclass' implementation of this method. this should
  // copy the output requested region to the input requested region
  Superclass::GenerateInputRequestedRegion();

  // This filter needs all of the input
  typename StructureTensorImageFilter<TInputImage,TOutputImage>
    ::InputImagePointer image = 
      const_cast<InputImageType *>( this->GetInput() );
  image->SetRequestedRegion( this->GetInput()->GetLargestPossibleRegion() );
}

template <typename TInputImage, typename TOutputImage, typename TOutputImage2 >
void
StructureTensorImageFilter<TInputImage,TOutputImage,TOutputImage2>
::EnlargeOutputRequestedRegion(DataObject *output)
{
  TOutputImage *out = dynamic_cast<TOutputImage*>(output);

  if ( out )
    {
    out->SetRequestedRegion( out->GetLargestPossibleRegion() );
    }
}

/**
 *   Make Ouput
 * \todo Verify that MakeOutput is createing the right type of objects
 *  this could be the cause of the reinterpret_cast bug in this class
 */
template <typename TInputImage, typename TOutputImage, typename TOutputImage2 >
DataObject::Pointer
StructureTensorImageFilter<TInputImage,TOutputImage,TOutputImage2>
::MakeOutput(unsigned int idx)
{
  DataObject::Pointer output;
  switch( idx )
    {
    case 0:
      output = (TOutputImage::New()).GetPointer();
      break;
    case 1:
      output = (TOutputImage2::New()).GetPointer();
      break;
    }
  return output.GetPointer();
}

template <typename TInputImage, typename TOutputImage, typename TOutputImage2 >
typename TOutputImage2::Pointer 
StructureTensorImageFilter<TInputImage,TOutputImage,TOutputImage2>
::GetSmoothGradient()
{
  return dynamic_cast< GradientImageType * >(
    this->ProcessObject::GetOutput(1) );
}

template <typename TInputImage, typename TOutputImage, typename TOutputImage2 >
typename TOutputImage::Pointer 
StructureTensorImageFilter<TInputImage,TOutputImage,TOutputImage2>
::GetTensor()
{
  return dynamic_cast< TensorImageType * >(
    this->ProcessObject::GetOutput(0));
}

/**
 * Compute filter for Gaussian kernel
 */
template <typename TInputImage, typename TOutputImage, typename TOutputImage2 >
void
StructureTensorImageFilter<TInputImage,TOutputImage,TOutputImage2>
::GenerateData(void)
{

  const typename TInputImage::ConstPointer inputImage( this->GetInput() );
  m_Tensor = this->GetTensor();
  m_Tensor->SetRegions(inputImage->GetRequestedRegion());
  m_Tensor->Allocate();
  m_Gradient = this->GetSmoothGradient();
  m_Gradient->SetRegions(inputImage->GetRequestedRegion());
  m_Gradient->Allocate();

  // Create a process accumulator for tracking the progress of this
  // minipipeline
  ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
  progress->SetMiniPipelineFilter(this);

  unsigned int imageDimensionMinus1 = static_cast<int>(ImageDimension)-1;
  // Compute the contribution of each filter to the total progress.
  const double weight = 1.0 / ( ImageDimension * ImageDimension );

  for( unsigned int i = 0; i<ImageDimension; i++ )
    {
    progress->RegisterInternalFilter( m_SmoothingFilters[i], weight );
    }
  progress->ResetProgress();

  // Step 1: Calculate gradient images using 1st order derivated gaussians

  for( size_t i=0; i<ImageDimension; ++i)
    {
    m_SmoothingFilters[ i ]->SetInput(inputImage);
    m_SmoothingFilters[ i ]->Update();
    }

  // Step 2: Tensor (outer, diadic) proouct of the gradient

  ImageRegionIterator< TensorImageType > ottensor( m_Tensor, 
    m_Tensor->GetRequestedRegion() );
  ImageRegionIterator< GradientImageType > otgradient( m_Gradient,
    m_Gradient->GetRequestedRegion() );
  std::vector<ImageRegionConstIterator<RealImageType> > 
    itgradient(ImageDimension);
  
  for( size_t i = 0;  i< ImageDimension; ++i)
    {
    itgradient[i] = ImageRegionConstIterator<RealImageType>(
    m_SmoothingFilters[i]->GetOutput(),
    m_SmoothingFilters[i]->GetOutput()->GetRequestedRegion());
    itgradient[i].GoToBegin();
    }

  ottensor.GoToBegin(); 
  otgradient.GoToBegin();
  GradientPixelType gradpx;
  TensorPixelType tenspx;
  while( !ottensor.IsAtEnd() )
    {
    unsigned int count = 0;
    for ( unsigned int j = 0; j < ImageDimension; ++j)
      {
      for (int k = j; k < ImageDimension; ++k)
        {
        tenspx[count++] = -itgradient[j].Get() * -itgradient[k].Get();
        }
      gradpx[j] = -itgradient[j].Get();
      ++itgradient[j];
      }
    ottensor.Set(tenspx);
    otgradient.Set(gradpx);
    ++ottensor;
    ++otgradient;
    }

    // Step 3: Integration filter (gaussian) in the tensor
    // Use RecursiveGaussian, DiscreteGaussianDerivative cannot handle tensors 
    typedef RecursiveGaussianImageFilter< TensorImageType, TensorImageType>
      TensorSmoothFilter;
    TensorSmoothFilter::Pointer IntegrationFilter = TensorSmoothFilter::New();
    IntegrationFilter->SetInput( m_Tensor );
    IntegrationFilter->SetOrder( TensorSmoothFilter::ZeroOrder );
    IntegrationFilter->SetSigma( m_rho );
    IntegrationFilter->Update();
    m_Tensor = IntegrationFilter->GetOutput();
}

template <typename TInputImage, typename TOutputImage, typename TOutputImage2>
void
StructureTensorImageFilter
  <TInputImage, TOutputImage, TOutputImage2>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  os << "NormalizeAcrossScale: " << m_NormalizeAcrossScale << "Sigma: " <<
  m_sigma << "Rho: " << m_rho << std::endl;
}

} // end namespace itk

#endif
