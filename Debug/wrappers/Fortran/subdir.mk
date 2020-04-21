################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../wrappers/Fortran/mui_3df.cpp 

OBJS += \
./wrappers/Fortran/mui_3df.o 

CPP_DEPS += \
./wrappers/Fortran/mui_3df.d 


# Each subdirectory must supply rules for building sources it contributes
wrappers/Fortran/%.o: ../wrappers/Fortran/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


