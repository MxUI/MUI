################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../wrappers/Python/mui4py/mui4py.cpp 

OBJS += \
./wrappers/Python/mui4py/mui4py.o 

CPP_DEPS += \
./wrappers/Python/mui4py/mui4py.d 


# Each subdirectory must supply rules for building sources it contributes
wrappers/Python/mui4py/%.o: ../wrappers/Python/mui4py/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


