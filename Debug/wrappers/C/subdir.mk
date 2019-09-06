################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../wrappers/C/mui_3d.cpp 

C_SRCS += \
../wrappers/C/unit_test.c 

OBJS += \
./wrappers/C/mui_3d.o \
./wrappers/C/unit_test.o 

CPP_DEPS += \
./wrappers/C/mui_3d.d 

C_DEPS += \
./wrappers/C/unit_test.d 


# Each subdirectory must supply rules for building sources it contributes
wrappers/C/%.o: ../wrappers/C/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

wrappers/C/%.o: ../wrappers/C/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


