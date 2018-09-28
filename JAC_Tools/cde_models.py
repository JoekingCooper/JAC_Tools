import numpy as np

class CDE_Modelling:
    """
    
        Class to model CDE lightcurves.
        
        This model uses a pixel based scheme to calculate each lightcurve and can be very slow to iterate over every pixel
        or have very poor resolution on the lightcurve.
        
        I recommend using the slice_lightcurve function, this is a good comprimise between the fastest model (which is the intensity of a single pixel)
        vs the slowest model which calculates the combined intensity from the star. This method can also give good estimates of the planetary radius 
        ratio when using a low resolution model, the resolution is set by the radius of the star and of the dust tail, so a radius < 40 pixels is 
        recommended.        
    
    """
    def __init__(self,prr=0.1,dtr=45,dc=-45,ipr=0.0,ii=160.0,ldca=0.3,ldcb=0.1,t0=0.0,sma=15,per=0.85):
        """
            Initialize the parameters of the class, mostly just stellar parameters. 
            
            Any of these can be changed, but I do not reccomend changing the resolution width and resolution height manually.
        """
        self.planetary_radius=prr
        self.dust_tail_ratio=dtr
        self.decay_constant=dc
        self.initial_intensity=ii
        self.impact_parameter=ipr
        self.limb_darkening_coeff_a=ldca
        self.limb_darkening_coeff_b=ldcb
        self.epoch_shift=t0
        self.period=per#in days
        self.semi_major=sma#needs to be in stellar radii
        
        self.dust_tail_length=self.planetary_radius*self.dust_tail_ratio
        self.star_horizontal_position=0.0
        self.star_vertical_position=0.0
        self.planet_vertical_position=self.star_vertical_position+self.impact_parameter
        self.orbital_speed=2*np.pi*self.semi_major/self.period
    def update_parameters(self):
        """
            If any parameters have been updated then please run this to update the other dependancies.
        """
        self.dust_tail_length=self.planetary_radius*self.dust_tail_ratio
        self.star_horizontal_position=0.0
        self.star_vertical_position=0.0
        self.planet_vertical_position=self.star_vertical_position+self.impact_parameter
        self.orbital_speed=2*np.pi*self.semi_major/self.period
    def I_star_prior(self,distance_ratio_from_stellar_centre):
        """
            Defining the boundry of the disc.
        """
        if distance_ratio_from_stellar_centre<1:
            return 1.0
        return 0.0
    def I_tail_prior(self,current_planet_horizontal_position,tail_lower_bound,tail_upper_bound,pixel_position_ver,pixel_position_hor):
        """
            Defining the boundry of the dust tail.
        """
        if current_planet_horizontal_position>pixel_position_hor>=current_planet_horizontal_position-self.dust_tail_length and tail_lower_bound<pixel_position_ver<tail_upper_bound:
            return 1.0
        return 0.0
    def I_core_prior(self,current_planet_horizontal_position,pixel_position_hor,pixel_position_ver):
        """
            Defining the boundry of the planetary core.
        """
        r_pl_temp=np.sqrt((self.planet_vertical_position-pixel_position_ver)**2+(current_planet_horizontal_position-pixel_position_hor)**2)
        if  r_pl_temp<self.planetary_radius and current_planet_horizontal_position<=pixel_position_ver<=current_planet_horizontal_position+self.planetary_radius:
            return 1.0
        return 0.0
    def I_star_antiprior(self,current_planet_horizontal_position,pixel_position_hor,pixel_position_ver,tail_lower_bound,tail_upper_bound):
        """
            Defines when the tail and core are present so that the stellar flux isn't added on top.
            
            There could be a more efficient way of doing this.
        """
        r_pl_temp=np.sqrt((self.planet_vertical_position-pixel_position_hor)**2+(current_planet_horizontal_position-pixel_position_ver)**2)
        #print(r_pl_comp_in,tail_l_in,x_mov,y_pos,i_in,y1_in,y2_in,m_in)
        if current_planet_horizontal_position>pixel_position_hor>=current_planet_horizontal_position-self.dust_tail_length and tail_lower_bound<pixel_position_ver<tail_upper_bound:
            return 0.0
        elif r_pl_temp<self.planetary_radius and current_planet_horizontal_position<=pixel_position_hor<=current_planet_horizontal_position+self.planetary_radius:
            return 0.0
        return 1.0
    def I_star_FUNC(self,distance_ratio_from_stellar_centre):
        """
            Equation defining the intensity of the star at a particular position from the centre of the stellar disc.
        """
        if distance_ratio_from_stellar_centre<1:
            mu_in=np.sqrt(1-distance_ratio_from_stellar_centre**2) 
            I_ratio=1-self.limb_darkening_coeff_a*(1-mu_in)-self.limb_darkening_coeff_b*(1-mu_in)**2
            #print(I_ratio,I0,r_in,r_comp_in,coeff_a,coeff_b,x_in,mu_in)
            I=self.initial_intensity*I_ratio
            return I
        return 0.0
    def I_core_FUNC(self,stellar_pixel_intensity,distance_from_planetary_centre):
        """ 
            Equation defining the intensity at distance from the planetary core.
        """
        trans_coeff=stellar_pixel_intensity/(np.log(self.planetary_radius+1))
        I=trans_coeff*np.log(distance_from_planetary_centre+1)
        if np.isnan(I):
            I=0.0
        return I
    def I_tail_FUNC(self,stellar_pixel_intensity,current_planet_horizontal_position,tail_lower_bound,pixel_position_hor,pixel_position_ver):
        """
            Equation defining the intensity of the tail at a distance from the tail centre and from the planetary core.
        """
        distance_from_tail_centre=abs(self.planet_vertical_position-pixel_position_ver)
        max_distance_from_tail_centre=abs(self.planet_vertical_position-tail_lower_bound)
        I_tail=(np.exp(-pixel_position_hor/self.decay_constant)-np.exp(-current_planet_horizontal_position/self.decay_constant))/(np.exp(-(current_planet_horizontal_position-self.dust_tail_length)/self.decay_constant)-np.exp(-(current_planet_horizontal_position/self.decay_constant)))
        trans_coeff=1/(np.log(max_distance_from_tail_centre+1))
        I_vert=trans_coeff*np.log(distance_from_tail_centre+1)
        I_combine=np.sqrt((I_vert)**2+I_tail**2)
        if I_combine>1.0:
            I_combine=1.0
        I=stellar_pixel_intensity*I_combine
        if np.isnan(I):
            I=0.0
        return I
    def uniformLimbDarkening_lightcurve_intensity_correction(self,quick_lightcurve):
        delta_flux=(max(quick_lightcurve)-min(quick_lightcurve))
        
        stellar_intensity=delta_flux/(self.planetary_radius**2)
        return quick_lightcurve+stellar_intensity-delta_flux
    def quick_lightcurve(self,t):
        """
            Lightcurve of the CDE for a single pixel at the centre of the star. This model gives no real information on the planetary radius.
        """
        quick_cde_lightcurve=[]
        for current_planet_horizontal_position in (t-self.epoch_shift)*self.orbital_speed:
            pixel_position_hor=self.star_horizontal_position 
            pixel_position_ver=self.star_vertical_position
            distance_from_planetary_centre=np.sqrt((self.star_vertical_position-pixel_position_ver)**2+(current_planet_horizontal_position-pixel_position_hor)**2)
            distance_ratio_from_stellar_centre=np.sqrt((self.star_vertical_position-pixel_position_ver)**2+(self.star_horizontal_position-pixel_position_hor)**2)
            tail_lower_bound=self.star_vertical_position-self.planetary_radius
            tail_upper_bound=self.star_vertical_position+self.planetary_radius
            stellar_pixel_intensity=self.I_star_prior(distance_ratio_from_stellar_centre)*self.I_star_FUNC(distance_ratio_from_stellar_centre)
            pixel_intensity=stellar_pixel_intensity*self.I_star_antiprior(current_planet_horizontal_position,pixel_position_hor,pixel_position_ver,tail_lower_bound,tail_upper_bound)+self.I_core_prior(current_planet_horizontal_position,pixel_position_hor,pixel_position_ver)*self.I_core_FUNC(stellar_pixel_intensity,distance_from_planetary_centre)+self.I_tail_FUNC(stellar_pixel_intensity,current_planet_horizontal_position,tail_lower_bound,pixel_position_hor,pixel_position_ver)*self.I_tail_FUNC(stellar_pixel_intensity,current_planet_horizontal_position,tail_lower_bound,pixel_position_hor,pixel_position_ver)
            quick_cde_lightcurve.append(pixel_intensity)
        quick_cde_lightcurve=self.uniformLimbDarkening_lightcurve_intensity_correction(quick_cde_lightcurve)
        quick_cde_lightcurve=np.array(quick_cde_lightcurve)/max(quick_cde_lightcurve)
        return quick_cde_lightcurve
    def slice_lightcurve(self,t,frame_steps=30):
        """
            Lightcurve of the CDE for a horizontal slice across the centre of the star. This model gives a good approximation of stellar/planetary parameters.
        """
        cde_lightcurve=[]
        pixel_position_hor_range=np.linspace(self.star_horizontal_position-1,self.star_horizontal_position+1,frame_steps)
        for current_planet_horizontal_position in (t-self.epoch_shift)*self.orbital_speed:
            frame_intensity=0
            for pixel_position_hor in pixel_position_hor_range:
                pixel_position_ver=self.star_vertical_position
                distance_from_planetary_centre=np.sqrt((self.star_vertical_position-pixel_position_ver)**2+(current_planet_horizontal_position-pixel_position_hor)**2)
                distance_ratio_from_stellar_centre=np.sqrt((self.star_vertical_position-pixel_position_ver)**2+(self.star_horizontal_position-pixel_position_hor)**2)
                tail_lower_bound=self.star_vertical_position-self.planetary_radius
                tail_upper_bound=self.star_vertical_position+self.planetary_radius
                stellar_pixel_intensity=self.I_star_prior(distance_ratio_from_stellar_centre)*self.I_star_FUNC(distance_ratio_from_stellar_centre)
                frame_intensity+=stellar_pixel_intensity*self.I_star_antiprior(current_planet_horizontal_position,pixel_position_hor,pixel_position_ver,tail_lower_bound,tail_upper_bound)+self.I_core_prior(current_planet_horizontal_position,pixel_position_hor,pixel_position_ver)*self.I_core_FUNC(stellar_pixel_intensity,distance_from_planetary_centre)+self.I_tail_FUNC(stellar_pixel_intensity,current_planet_horizontal_position,tail_lower_bound,pixel_position_hor,pixel_position_ver)*self.I_tail_FUNC(stellar_pixel_intensity,current_planet_horizontal_position,tail_lower_bound,pixel_position_hor,pixel_position_ver)
            cde_lightcurve.append(frame_intensity)
        cde_lightcurve=np.array(cde_lightcurve)/max(cde_lightcurve)
        return cde_lightcurve
    def full_lightcurve(self,t,frame_steps=50):
        """
            Lightcurve of the CDE for the entire frame. This model gives a very good approximation of stellar/planetary parameters.
        """
        cde_lightcurve=[]
        pixel_position_hor_range=np.linspace(self.star_horizontal_position-1,self.star_horizontal_position+1,frame_steps)
        pixel_position_ver_range=np.linspace(self.star_vertical_position-1,self.star_vertical_position+1,frame_steps)
        for current_planet_horizontal_position in (t-self.epoch_shift)*self.orbital_speed:
            frame_intensity=0
            for pixel_position_hor in pixel_position_hor_range:
                for pixel_position_ver in pixel_position_ver_range:
                    distance_from_planetary_centre=np.sqrt((self.star_vertical_position-pixel_position_ver)**2+(current_planet_horizontal_position-pixel_position_hor)**2)
                    distance_ratio_from_stellar_centre=np.sqrt((self.star_vertical_position-pixel_position_ver)**2+(self.star_horizontal_position-pixel_position_hor)**2)
                    tail_lower_bound=self.star_vertical_position-self.planetary_radius
                    tail_upper_bound=self.star_vertical_position+self.planetary_radius
                    stellar_pixel_intensity=self.I_star_prior(distance_ratio_from_stellar_centre)*self.I_star_FUNC(distance_ratio_from_stellar_centre)
                    pixel_intensity=stellar_pixel_intensity*self.I_star_antiprior(current_planet_horizontal_position,pixel_position_hor,pixel_position_ver,tail_lower_bound,tail_upper_bound)+self.I_core_prior(current_planet_horizontal_position,pixel_position_hor,pixel_position_ver)*self.I_core_FUNC(stellar_pixel_intensity,distance_from_planetary_centre)+self.I_tail_FUNC(stellar_pixel_intensity,current_planet_horizontal_position,tail_lower_bound,pixel_position_hor,pixel_position_ver)*self.I_tail_FUNC(stellar_pixel_intensity,current_planet_horizontal_position,tail_lower_bound,pixel_position_hor,pixel_position_ver)
                    frame_intensity+=pixel_intensity
            cde_lightcurve.append(frame_intensity)
        cde_lightcurve=np.array(cde_lightcurve)/max(cde_lightcurve)
        return cde_lightcurve


"""
#This is just an example of what you can do
t=np.linspace(-0.1,0.1,100)
CDEM=CDE_Modelling(dtr=10,prr=0.6,t0=-0.01,sma=200)
y=CDEM.quick_lightcurve(t)
#y=CDEM.full_lightcurve(frame_steps=30)
#y=CDEM.slice_lightcurve(t)
import matplotlib.pyplot as plt
plt.figure()
plt.plot(t,y)
#print(CDEM.resolution_width,CDEM.star_horizontal_position)
plt.show()
"""