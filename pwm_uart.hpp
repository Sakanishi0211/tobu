#ifndef PWM_UART_HPP
#define PWM_UART_HPP
#include <stdio.h>
#include "hardware/irq.h"
#include "hardware/uart.h"
#include "pico/stdlib.h"
#include "hardware/pwm.h"
#define dutymin 1330
#define dutymax 1520
#define UART_ID uart0
#define BAUD_RATE 100000      //115200
#define DATA_BITS 8
#define STOP_BITS 2             //1
#define PARITY    UART_PARITY_EVEN
#define DUTYMIN 1330
#define DUTYMAX 1520
#define CH1MAX 1680
#define CH1MIN 368
#define CH2MAX 1680
#define CH2MIN 368
#define CH3MAX 1680
#define CH3MIN 368
#define CH4MAX 1680
#define CH4MIN 368
#define CH5MAX 1904
#define CH5MIN 144
#define CH6MAX 1904
#define CH6MIN 144
#define UART_TX_PIN 0
#define UART_RX_PIN 1
extern float Data1,Data2,Data3,Data4,Data5,Data6;
extern uint slice_num[2];
uint8_t serial_settei(void);
uint8_t pwm_settei(void);

#endif
