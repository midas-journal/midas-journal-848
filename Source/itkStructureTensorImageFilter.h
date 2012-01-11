
#ifndef __itkStructureTensorImageFilter_h
#define __itkStructureTensorImageFilter_h

#include "itkRecursiveGaussianImageFilter.h"
#include "itkDiscreteGaussianDerivativeImageFilter.h"
#include "itkNthElementImageAdaptor.h"
#include "itkImage.h"
#include "itkSymmetricSecondRankTensor.h"
#include "itkPixelTraits.h"
#include "itkProgressAccumulator.h"

namespace itk
{

/** \class StructureTensorImageFilter
 *
 * \brief Computes the structure tensor of a multidimensional image 
 * 
 * 
 * \ingroup GradientFilters
 * \ingroup Singlethreaded
 */
// NOTE that the ITK_TYPENAME macro has to be used here in lieu 
// of "typename" because VC++ doesn't like the typename keyword 
// on the defaults of template parameters

template <typename TInputImage,
          typename TOutputImage= Image< SymmetricSecondRankTensor< 
            ITK_TYPENAME NumericTraits< ITK_TYPENAME TInputImage::PixelType>
            ::RealType,
            ::itk::GetImageDimension<TInputImage>::ImageDimension >,
            ::itk::GetImageDimension<TInputImage>::ImageDimension >,
          typename TOutputImage2 = Image< CovariantVector< 
            ITK_TYPENAME NumericTraits< ITK_TYPENAME TInputImage::PixelType>
            ::RealType, 
            ::itk::GetImageDimension<TInputImage>::ImageDimension >, 
            ::itk::GetImageDimension<TInputImage>::ImageDimension > >
class ITK_EXPORT StructureTensorImageFilter:
    public ImageToImageFilter<TInputImage,TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef StructureTensorImageFilter                   Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage> Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;
  
  /** Pixel Type of the input image */
  typedef TInputImage                                    InputImageType;
  typedef typename TInputImage::PixelType                PixelType;
  typedef typename NumericTraits<PixelType>::RealType    RealType;


  /** Image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Define the image type for internal computations 
      RealType is usually 'double' in NumericTraits. 
      Here we prefer float in order to save memory.  */

  typedef float                              InternalRealType;
  typedef Image<InternalRealType, 
    itkGetStaticConstMacro(ImageDimension) > RealImageType;

  /**  Smoothing filter type */
  
  typedef DiscreteGaussianDerivativeImageFilter<RealImageType, RealImageType>
    GaussianFilterType;

  /**  Pointer to a gaussian filter.  */
  typedef typename GaussianFilterType::Pointer    GaussianFilterPointer;

  /** Type of the output Image */
  typedef TOutputImage                                     TensorImageType;
  typedef TOutputImage2                                    GradientImageType;
  typedef typename          TensorImageType::PixelType     TensorPixelType;
  typedef typename          GradientImageType::PixelType   GradientPixelType;
  typedef typename          TensorImageType::Pointer       TensorImagePointer;
  typedef typename          GradientImageType::Pointer     GradientImagePointer;
  typedef typename PixelTraits<TensorPixelType>::ValueType OutputComponentType;
  typedef itk::Image< CovariantVector< InternalRealType, ImageDimension >, 
    ImageDimension > VectorImageType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Set Sigma value. Sigma is measured in the units of image spacing.  */
  void SetSigma( RealType sigma );

  /** Set Rho value. Sigma is measured in the units of image spacing.  */
  void SetRho( RealType rho );

  /** Define which normalization factor will be used for the Gaussian */
  void SetNormalizeAcrossScale( bool normalizeInScaleSpace );
  itkGetMacro( NormalizeAcrossScale, bool );

  /** StructureTensorImageFilter needs all of the input to produce an
   * output. Therefore, StructureTensorImageFilter needs to provide
   * an implementation for GenerateInputRequestedRegion in order to inform
   * the pipeline execution model.
   * \sa ImageToImageFilter::GenerateInputRequestedRegion() */
  virtual void GenerateInputRequestedRegion() 
    throw(InvalidRequestedRegionError);

  /**  Create the Output */
  DataObject::Pointer MakeOutput(unsigned int idx);

  GradientImagePointer GetSmoothGradient(void);
  TensorImagePointer GetTensor(void);

protected:
  StructureTensorImageFilter();
  virtual ~StructureTensorImageFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;
  
  /** Generate Data */
  void GenerateData( void );

  // Override since the filter produces the entire dataset
  void EnlargeOutputRequestedRegion(DataObject *output);

private:
  StructureTensorImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  std::vector<GaussianFilterPointer>         m_SmoothingFilters;
  TensorImagePointer                         m_Tensor;
  GradientImagePointer                       m_Gradient;

  /** Normalize the image across scale space */
  bool     m_NormalizeAcrossScale; 
  RealType m_sigma;
  RealType m_rho;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkStructureTensorImageFilter.txx"
#endif

#endif
